import rasterio
import rasterio.features
import rioxarray
import geopandas
import pandas
import numpy as np
import scipy.ndimage
import scipy.interpolate
import whitebox
import pysheds
import pysheds.grid
import pysheds.view
import skimage.feature
import shapely
import tempfile
import logging
import os
from tqdm import tqdm
import argparse


def extract_ice_divides(
    input_outlines_path,
    input_dem_path,
    output_path,
    fill_holes=False,
    dem_conditioned=False,
    dem_conditioning_method="breach",
    dem_smoothing=7,
    pixel_size=10.0,
    buffer1_size=160.0,
    buffer2_size=60.0,
    gutter_size=1,
    l_gutter=100.0,
    routing="dinf",
    sortby="accumulation",
    a=5.72,
    b=0.5,
    c=1400,
    percentile=0.80,
    skip_missing_artifacts=False,
    a_sliver_1=1e5, 
    a_sliver_2=2e5, 
    c_sliver=0.5,
):
    """
    TODO: __docstring__
    """
    if not routing in {"d8", "mfd", "dinf"}:
        raise ValueError()
    if not sortby in {"accumulation", "elevation"}:
        raise ValueError()
    buffer1_size = int(buffer1_size / pixel_size + 0.5)
    buffer2_size = int(buffer2_size / pixel_size + 0.5)
    a, c = a / pixel_size, c / pixel_size
    shed_kwargs = get_shed_kwargs(routing)
    maxima_search_radius = 4 * gutter_size + 1
    accumulation_threshold = (buffer2_size / a)**(1/b)
    
    outlines = get_conditioned_outlines(input_outlines_path, fill_holes)
    output = np.zeros(outlines.shape, dtype=np.uint32)
    _, breached_dem_path = get_conditioned_dem(input_dem_path, dem_conditioned, dem_conditioning_method, dem_smoothing)
    complexes, n_complexes = scipy.ndimage.label(outlines, structure=np.ones((3, 3)))
    
    for complex_idx in range(1, n_complexes + 1):
        # if complex_idx not in [1, 269]:
        #     continue
        print(f"Processing ice complex {complex_idx}/{n_complexes}")
        # extract ice complex and make buffers
        complex = (complexes == complex_idx)
        complex, (complex_min_row, complex_max_row, complex_min_col, complex_max_col) = \
            trim_mask(complex, buffer_size=buffer1_size, fillvalue=0)
        buffer1 = scipy.ndimage.binary_dilation(complex, iterations=buffer1_size)
        buffer2 = scipy.ndimage.binary_dilation(complex, iterations=buffer2_size)
        gutter = scipy.ndimage.binary_dilation(complex, iterations=buffer2_size + gutter_size)
        gutter[buffer2 == 1] = 0
        
        # extract dem for ice complex
        grid = pysheds.grid.Grid.from_raster(breached_dem_path)
        complex_dem = grid.read_raster(breached_dem_path)
        complex_dem = complex_dem[complex_min_row:complex_max_row + 1, complex_min_col:complex_max_col + 1]
        grid.clip_to(complex_dem)
        complex_dem[buffer1 == 0] = 9999
        complex_dem[gutter == 1] -= l_gutter
        
        # calculate flow accumulation
        fdir = grid.flowdir(complex_dem, **shed_kwargs)
        acc = grid.accumulation(fdir, **shed_kwargs)
        # extract potential pour points
        local_maxima = extract_local_maxima(acc, gutter == 1, maxima_search_radius)
        potential_pour_points = [(row, col, acc[row, col], complex_dem[row, col]) for row, col in local_maxima]
        potential_pour_points = clean_potential_pour_point(
            potential_pour_points, accumulation_threshold, sortby, percentile
        )
        
        # get flowsheds and pour points, merge flowsheds
        flowsheds, pour_points = get_flowsheds_and_pour_points(complex, potential_pour_points, grid, fdir, shed_kwargs)
        intersections = find_intersections(flowsheds, pour_points, complex, buffer2_size + gutter_size, a, b, c)
        complex_glaciers = np.zeros(complex.shape, dtype=np.uint32)
        glacier_idx = 1
        for intersection in intersections:
            for flowshed_idx in intersection:
                complex_glaciers[flowsheds == flowshed_idx] = glacier_idx
            glacier_idx += 1
        missed = (complex == 1) & (complex_glaciers == 0)
        if not skip_missing_artifacts and missed.any():
            logging.warning("Some pixels of ice complex were missed. They will be forcely written into output. Thus, reliable ice divides are not guaranteed. However, it could be isolated small polygon. It may be useful to consider choosing different set of algorithm parameters.")
            missed, n_missed = scipy.ndimage.label(missed, structure=np.ones((3, 3)))
            for missed_idx in range(1, n_missed + 1):
                complex_glaciers[missed == missed_idx] = glacier_idx
                glacier_idx += 1
        # map complex_glaciers back to global scene
        output[complex_min_row:complex_max_row + 1, complex_min_col:complex_max_col + 1] += complex_glaciers
    
    # polygonize, filter out slivers and save
    glaciers = polygonize(output, output_path, input_outlines_path)
    glaciers = filter_out_slivers(glaciers, a_sliver_1, a_sliver_2, c_sliver)
    glaciers = glaciers.set_crs(outlines.rio.crs)
    glaciers.to_file(output_path)
    # clean temporary
    if breached_dem_path != input_dem_path:
        os.remove(breached_dem_path)


def get_shed_kwargs(routing):
    if routing == "d8":
        dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
        shed_kwargs = {"dirmap": dirmap, "routing": "d8"}
    elif routing == "mfd":
        shed_kwargs = {"routing": "mfd"}
    elif routing == "dinf":
        shed_kwargs = {"routing": "dinf"}
    else:
        raise NotImplementedError()
    return shed_kwargs


def get_conditioned_outlines(outlines_path, fill_holes):
    conditioned_outlines = rioxarray.open_rasterio(outlines_path)[0]
    if not fill_holes:
        return conditioned_outlines
    conditioned_outlines = scipy.ndimage.binary_fill_holes(conditioned_outlines)
    return conditioned_outlines
    

def get_conditioned_dem(dem_path, conditioned, method, smoothing):
    if conditioned:
        return rioxarray.open_rasterio(dem_path)[0], dem_path
    
    wbt = whitebox.WhiteboxTools()
    
    conditioned_dem = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
    conditioned_dem_path = conditioned_dem.name
    
    if smoothing:
        wbt.feature_preserving_smoothing(dem_path, conditioned_dem_path, filter=smoothing)
        dem_path = conditioned_dem_path
    
    if method == "fill":
        wbt.fill_single_cell_pits(dem_path, conditioned_dem_path)
        wbt.fill_depressions(conditioned_dem_path, conditioned_dem_path)
    elif method == "breach":
        wbt.breach_single_cell_pits(dem_path, conditioned_dem_path)
        wbt.breach_depressions(conditioned_dem_path, conditioned_dem_path)
    else:
        raise NotImplementedError()
        
    conditioned_dem.close()
    conditioned_dem = rioxarray.open_rasterio(conditioned_dem_path)[0]
    return conditioned_dem, conditioned_dem_path


def trim_mask(mask, buffer_size=0, fillvalue=0):
    height, width = mask.shape
    rows, cols = np.where(mask != fillvalue)
    min_row, max_row = np.min(rows), np.max(rows)
    min_col, max_col = np.min(cols), np.max(cols)
    min_row, max_row = min_row - buffer_size, max_row + buffer_size
    min_col, max_col = min_col - buffer_size, max_col + buffer_size
    min_row, min_col = max(min_row, 0), max(min_col, 0)
    max_row, max_col = min(max_row, height - 1), min(max_col, width - 1)
    trimmed = mask[min_row:max_row + 1, min_col:max_col + 1]
    return trimmed, (min_row, max_row, min_col, max_col)
    
    
def extract_local_maxima(arr, mask, search_radius):
    # dem conditioning may alter the shape...
    if arr.shape != mask.shape:
        new_height = min(arr.shape[0], mask.shape[0])
        new_width = min(arr.shape[1], mask.shape[1])
        arr = arr[:new_height, :new_width]
        mask = mask[:new_height, :new_width]
    footprint = np.ones((search_radius, search_radius))
    local_maxima_arr = mask & \
        (scipy.ndimage.maximum_filter(arr, footprint=footprint, mode="constant") == arr)
    local_maxima_rows, local_maxima_cols = np.where(local_maxima_arr == 1)
    local_maxima = list(zip(local_maxima_rows, local_maxima_cols))
    return local_maxima
    
    
def clean_potential_pour_point(potential_pour_points, accumulation_threshold, sortby, percentile):
    potential_pour_points = [_ for _ in potential_pour_points if _[2] > accumulation_threshold]
    if sortby == "accumulation":
        key = lambda x: x[2]
    elif sortby == "elevation":
        key = lambda x: -x[3]
    else:
        raise NotImplementedError()
    potential_pour_points.sort(key=key)
    n_potential_pour_points = len(potential_pour_points)
    potential_pour_points = potential_pour_points[int(percentile * n_potential_pour_points):]
    return potential_pour_points


def get_flowsheds_and_pour_points(complex, potential_pour_points, grid, fdir, shed_kwargs):
    # sometimes pysheds trims complex..., so trim both to match
    # if complex.shape != fdir.shape:
    #     new_height = min(complex.shape[0], fdir.shape[0])
    #     new_width = min(complex.shape[1], fdir.shape[1])
    #     complex = complex[:new_height, :new_width]
    #     fdir = fdir[:new_height, :new_width]

    height, width = complex.shape
    watersheds = np.zeros((height, width), dtype=np.uint32)
    pour_points = dict()
    watershed_idx = 1

    with tqdm(total=len(potential_pour_points), desc="Calculating flowsheds") as pbar:
        while potential_pour_points:
            pour_point = potential_pour_points.pop()
            row, col, _, _ = pour_point
            if watersheds[row, col]:
                pbar.update(1)
                continue
            watershed = grid.catchment(x=col, y=row, fdir=fdir, xytype="index", **shed_kwargs)

            pour_points_to_delete = set()
            for other_pour_point, other_watershed_idx in pour_points.items():
                other_row, other_col, _, _ = other_pour_point
                if watershed[other_row, other_col]:
                    watershed[watersheds == other_watershed_idx] = True
                    pour_points_to_delete.add(other_pour_point)

            pour_points[pour_point] = watershed_idx

            for pour_point_to_delete in pour_points_to_delete:
                del pour_points[pour_point_to_delete]

            # dem conditioning may alter the shape...
            if watersheds.shape != watershed.shape:
                curated_watershed = np.zeros_like(watersheds, dtype=bool)
                max_height, max_width = min(height, watershed.shape[0]), min(width, watershed.shape[1])
                curated_watershed[:max_height, :max_width] = watershed[:max_height, :max_width]
                watershed = curated_watershed
            watersheds[watershed] = watershed_idx
            watershed_idx += 1

            pbar.update(1)

    flowsheds = watersheds.copy()
    flowsheds[complex == 0] = 0
    
    # interpolate only with adjacent flowsheds, otherwise it produces artifacts
    gaps = (complex == 1) & (flowsheds == 0)
    gaps, n_gaps = scipy.ndimage.label(gaps)
    curated_flowsheds = flowsheds.copy()

    for gap_idx in tqdm(range(1, n_gaps + 1), desc="Interpolating not delineated"):
        gap = (gaps == gap_idx)
        # make a subset for faster computations
        gap_subset, (min_row, max_row, min_col, max_col) = \
            trim_mask(gap, buffer_size=1, fillvalue=0)
        flowsheds_subset = flowsheds[min_row:max_row + 1, min_col:max_col + 1]
        dilated_gap = scipy.ndimage.binary_dilation(gap_subset, iterations=1)
        flowshed_idxs = np.unique(flowsheds_subset[(dilated_gap == 1) & (flowsheds_subset != 0)])
        if not len(flowshed_idxs):
            continue
        adjacent_flowsheds = np.isin(flowsheds_subset, flowshed_idxs)
        rows, cols = np.where(adjacent_flowsheds)
        values = flowsheds_subset[rows, cols]
        out_rows, out_cols = np.where(gap_subset)
        interpolated_values = scipy.interpolate.griddata((cols, rows), values, (out_cols, out_rows), method="nearest")
        curated_flowsheds[out_rows + min_row, out_cols + min_col] = interpolated_values
    
    return curated_flowsheds, pour_points
    

def find_intersections(flowsheds, pour_points, complex, buffer_size, a, b, c):
    intersections = []

    for pour_point in tqdm(pour_points, desc="Merging flowsheds"):
        row, col, _, _ = pour_point
        rows, cols = np.meshgrid(
            np.arange(flowsheds.shape[0]), np.arange(flowsheds.shape[1]), indexing="ij"
        )
        distance = np.sqrt((rows - row)**2 + (cols - col)**2)
        p_glac = distance <= min(c, a * pour_point[2] ** b)
        p_glac = p_glac & complex
        if not np.sum(p_glac != 0):
            continue
        p_glac, n_components = scipy.ndimage.label(p_glac, structure=np.ones((3, 3)))
        
        # make subset to speedup computations
        p_glac_subset, (min_row, max_row, min_col, max_col) = \
            trim_mask(p_glac, buffer_size=buffer_size, fillvalue=0)

        flowshed_idx = pour_points[pour_point]
        shortest_distance = np.inf
        closest_component = None
        for component_idx in range(1, n_components + 1):
            component = (p_glac_subset == component_idx)
            distances = scipy.ndimage.distance_transform_edt(~component)
            distance = distances[row - min_row, col - min_col]
            if distance < shortest_distance:
                shortest_distance = distance
                closest_component = component_idx

        p_glac = (p_glac == closest_component)
        intersected_flowsheds = frozenset(np.unique(flowsheds[p_glac]))
        intersections.append(intersected_flowsheds)
        
    grouped_intersections = union_intersecting(intersections)
    return grouped_intersections
    
    
def union_intersecting(sets):
    grouped_intersections = []
    visited = set()
    for s in sets:
        if s in visited:
            continue
        union = set(s)
        stack = [s]
        while stack:
            current_set = stack.pop()
            visited.add(current_set)
            for other_set in sets:
                if other_set in visited:
                    continue
                if not current_set.intersection(other_set):
                    continue
                union.update(other_set)
                stack.append(other_set)
                visited.add(other_set)
        grouped_intersections.append(union)
    return grouped_intersections


def polygonize(arr, output_path, reference_path):
    with rasterio.open(reference_path, "r") as src:
        trans = src.transform
        crs = src.crs

    polygons = [
        {"properties": {"id": label}, "geometry": geom}
        for geom, label in rasterio.features.shapes(arr.astype(np.int32), transform=trans)
        if label != 0
    ]

    dataframe = geopandas.GeoDataFrame.from_features(polygons, crs=crs)
    dataframe.to_file(output_path)
    return dataframe

    
def filter_out_slivers(glaciers, a_sliver_1, a_sliver_2, c_sliver):
    glaciers["area"] = glaciers.geometry.area
    glaciers["perimeter"] = glaciers.geometry.length
    glaciers["compactness"] = 2 * np.sqrt(glaciers.area * np.pi) / glaciers.perimeter
    
    slivers = glaciers[(glaciers.area < a_sliver_1) | ((glaciers.area < a_sliver_2) & (glaciers.compactness < c_sliver))]
    non_slivers = glaciers[~glaciers.index.isin(slivers.index)]
    
    # merge adjacent slivers
    polygons_to_merge = []
    for sliver in slivers.itertuples():
        neighbor_slivers = slivers[slivers.geometry.touches(sliver.geometry)]
        if neighbor_slivers.empty:
            continue
        polygons_to_merge.append(frozenset(neighbor_slivers.index) | {sliver.Index})
    polygons_to_merge = union_intersecting(polygons_to_merge)
    
    for s in polygons_to_merge:
        base_sliver_idx = None

        for sliver_idx in s:
            if base_sliver_idx is None:
                base_sliver_idx = sliver_idx
                continue
            geometry = slivers.loc[base_sliver_idx].geometry
            geometry = shapely.ops.unary_union([geometry, slivers.loc[sliver_idx].geometry])
            slivers.loc[base_sliver_idx, "geometry"] = geometry
            slivers = slivers[slivers.index != sliver_idx]
            
    # merge non-slivers with slivers
    polygons_to_merge = {}
    for sliver in slivers.itertuples():
        neighbors = non_slivers[non_slivers.geometry.touches(sliver.geometry)]
        if neighbors.empty:
            # do not filter out slivers not adjacent to other polygons
            # non_slivers = non_slivers.append(slivers.loc[sliver.Index], ignore_index=False)
            non_slivers = pandas.concat([non_slivers, geopandas.GeoDataFrame([slivers.loc[sliver.Index]])])
            slivers = slivers[slivers.index != sliver.Index]
            continue
        neighbor_with_longest_boundary = max(
            neighbors.itertuples(), 
            key=lambda x: x.geometry.intersection(sliver.geometry).length
        )
        if neighbor_with_longest_boundary.Index in polygons_to_merge:
            polygons_to_merge[neighbor_with_longest_boundary.Index].append(sliver.Index)
        else:
            polygons_to_merge[neighbor_with_longest_boundary.Index] = [sliver.Index]

    cleaned_glaciers = non_slivers.copy()
    for neighbor_idx, sliver_idxs in polygons_to_merge.items():
        geometry = non_slivers.loc[neighbor_idx].geometry
        for sliver_idx in sliver_idxs:
            geometry = shapely.ops.unary_union([geometry, slivers.loc[sliver_idx].geometry])
        cleaned_glaciers.loc[neighbor_idx, "geometry"] = geometry
    
    cleaned_glaciers = cleaned_glaciers.drop(columns=["area", "perimeter", "compactness"])
    return cleaned_glaciers

    
if __name__ == "__main__":  
    parser = argparse.ArgumentParser()
    
    parser.add_argument("input_outlines_path", help="Path to input outlines .tif file")
    parser.add_argument("input_dem_path", help="Path to input DEM .tif file")
    parser.add_argument("output_path", help = "Path to output .shp file")
    
    parser.add_argument("--fill_holes", action="store_true", help="Fill holes in the outlines")
    
    parser.add_argument("--dem_conditioned", action="store_true", help="Skip DEM conditioning")
    parser.add_argument("--dem_conditioning_method", default="fill", choices=["fill", "breach"], help="DEM conditioning method")
    parser.add_argument("--dem_smoothing", type=int, default=7, help="DEM smoothing parameter")
    
    parser.add_argument("--pixel_size", type=float, default=10.0, help="Pixel size in CRS units")
    parser.add_argument("--buffer1_size", type=float, default=160.0, help="Outer buffer size in CRS units")
    parser.add_argument("--buffer2_size", type=float, default=60.0, help="Inner buffer size in CRS units")
    parser.add_argument("--gutter_size", type=int, default=1, help="Gutter size in pixels")
    parser.add_argument("--l_gutter", type=float, default=100.0, help="Gutter depth")
    
    parser.add_argument("--routing", default="dinf", choices=["d8", "mfd", "dinf"], help="Routing method for hydrological calculations")
    parser.add_argument("--sortby", default="accumulation", choices=["accumulation", "elevation"], help="Sorting criteria for pour points")
    parser.add_argument("--a", type=float, default=5.72, help="Radii calculation parameter a")
    parser.add_argument("--b", type=float, default=0.5, help="Radii calculation parameter b")
    parser.add_argument("--c", type=float, default=1400.0, help="Maximum radius to merge pour points")
    parser.add_argument("--percentile", type=float, default=0.80, help="Percentile to filter out pour point with low accumulation/elevation")
    
    parser.add_argument("--skip_missing_artifacts", action="store_true", help="Skip pixels that were missed because of poor calibration/bad data/software errors")
    
    parser.add_argument("--a_sliver_1", type=float, default=1e5, help="Area threshold for filtering out slivers")
    parser.add_argument("--a_sliver_2", type=float, default=2e5, help="Area threshold for filtering out prolongated slivers")
    parser.add_argument("--c_sliver", type=float, default=0.5, help="Compactness threshold for filtering out slivers")
    
    args = parser.parse_args()
    
    extract_ice_divides(
        os.path.abspath(args.input_outlines_path),
        os.path.abspath(args.input_dem_path),
        os.path.abspath(args.output_path),
        fill_holes=args.fill_holes,
        dem_conditioned=args.dem_conditioned,
        dem_conditioning_method=args.dem_conditioning_method,
        dem_smoothing=args.dem_smoothing,
        pixel_size=args.pixel_size,
        buffer1_size=args.buffer1_size,
        buffer2_size=args.buffer2_size,
        gutter_size=args.gutter_size,
        l_gutter=args.l_gutter,
        routing=args.routing,
        sortby=args.sortby,
        a=args.a,
        b=args.b,
        c=args.c,
        percentile=args.percentile,
        a_sliver_1=args.a_sliver_1, 
        a_sliver_2=args.a_sliver_2, 
        c_sliver=args.c_sliver,
    )

import geopandas
import pandas
import rasterio
import rasterio.features
import rioxarray
import whitebox
import tempfile
# import cv2
import skimage.segmentation
import numpy as np
import shapely
import argparse
import os


def copy_ice_divides(
    ref_path, 
    complexes_path, 
    output_path, 
    dem_path,
    dem_conditioned=False,
    dem_conditioning_method="fill",
    dem_smoothing=7,
    no_sliver_filtering=False,
    a_sliver_1=1e5, 
    a_sliver_2=2e5, 
    c_sliver=0.5,
):
    """
    TODO: __docstring__
    """
    ref = geopandas.read_file(ref_path)
    complexes = geopandas.read_file(complexes_path)
    dem_arr, conditioned_dem_path = get_conditioned_dem(dem_path, dem_conditioned, dem_conditioning_method, dem_smoothing)
    with rasterio.open(conditioned_dem_path, "r") as dem_src:
        dem_transform = dem_src.transform
        dem_meta = dem_src.meta.copy()
        
    ref_shapes = ((geom, idx + 1) for idx, geom in enumerate(ref.geometry))
    complexes_shapes = ((geom, 1) for geom in complexes.geometry)
    
    ref_arr = rasterio.features.rasterize(shapes=ref_shapes, fill=0, out_shape=dem_arr.shape, transform=dem_transform)
    complexes_arr = rasterio.features.rasterize(shapes=complexes_shapes, fill=0, out_shape=dem_arr.shape, transform=dem_transform)
    mask = (ref_arr > 0) | (complexes_arr > 0)
        
    output = skimage.segmentation.watershed(-dem_arr, ref_arr, mask=mask)
    output[complexes_arr == 0] = 0 # trim to the original prediction
    output[(complexes_arr > 0) & (output == 0)] = 1 # restore isolated patches filtered out after watersheds
        
    # polygonize, filter out slivers and save
    glaciers = polygonize(output, output_path, dem_path)
    if not no_sliver_filtering:  
        glaciers = filter_out_slivers(glaciers, a_sliver_1, a_sliver_2, c_sliver)
    glaciers = glaciers.set_crs(ref.crs)
    glaciers.to_file(output_path)
    # remove temporary files
    if conditioned_dem_path != dem_path:
        os.remove(conditioned_dem_path)    
    

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_path", help="Path to reference .shp file")
    parser.add_argument("complexes_path", help="Path to to-be-divided .shp file")
    parser.add_argument("output_path", help="Path to output .shp file")
    parser.add_argument("dem_path", help="Path to elevation .tif file")
    
    parser.add_argument("--dem_conditioned", action="store_true", help="Skip DEM conditioning")
    parser.add_argument("--dem_conditioning_method", default="fill", choices=["fill", "breach"], help="DEM conditioning method")
    parser.add_argument("--dem_smoothing", type=int, default=7, help="DEM smoothing parameter")
    
    parser.add_argument("--no_sliver_filtering", action="store_true", help="Skip sliver filtering")
    parser.add_argument("--a_sliver_1", type=float, default=1e5, help="Area threshold for filtering out slivers")
    parser.add_argument("--a_sliver_2", type=float, default=2e5, help="Area threshold for filtering out prolongated slivers")
    parser.add_argument("--c_sliver", type=float, default=0.5, help="Compactness threshold for filtering out slivers")
    
    args = parser.parse_args()

    copy_ice_divides(
        os.path.abspath(args.ref_path), 
        os.path.abspath(args.complexes_path), 
        os.path.abspath(args.output_path), 
        os.path.abspath(args.dem_path),
        dem_conditioned=args.dem_conditioned,
        dem_conditioning_method=args.dem_conditioning_method,
        dem_smoothing=args.dem_smoothing,
        no_sliver_filtering=args.no_sliver_filtering,
        a_sliver_1=args.a_sliver_1, 
        a_sliver_2=args.a_sliver_2, 
        c_sliver=args.c_sliver,
    )
    
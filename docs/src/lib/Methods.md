```@raw html
<!--
  Library reference page.

  Each @docs block here lists bindings that have a docstring in
  src/PhysicalMeshes.jl (or in any of the included files). Use the
  index at the bottom of the page to jump to a specific name; use the
  manual pages for design rationale and extended examples.

  Generated content lives in this file — edit the manual pages, not
  this one — except for blocks here that have to be hand-curated to
  stay useful.
-->
```

# Library Reference

The complete API surface of `PhysicalMeshes.jl`. The manual pages discuss the
design rationale and walk through extended examples; this page is the
canonical "what does this name refer to?" reference. Bindings are grouped by
where they live in the source tree.

For a flat index, see the bottom of this page.

## Geometry primitives

### Concrete shapes

```@docs
PhysicalMeshes.Line2D
PhysicalMeshes.Line
PhysicalMeshes.Plane
PhysicalMeshes.Polygon2D
PhysicalMeshes.Polygon3D
PhysicalMeshes.Triangle
PhysicalMeshes.Triangle2D
PhysicalMeshes.Tetrahedron
PhysicalMeshes.Cube2D
PhysicalMeshes.Cube
PhysicalMeshes.Sphere
PhysicalMeshes.Ray2D
PhysicalMeshes.Ray3D
```

### Lines, planes, polygons, triangles, cubes, tetrahedra

```@docs
PhysicalMeshes.normal(::PhysicalMeshes.Plane)
PhysicalMeshes.distance
PhysicalMeshes.coplanar
PhysicalMeshes.polygon_rect
PhysicalMeshes.polygon_regular
PhysicalMeshes.area(::PhysicalMeshes.Polygon2D)
PhysicalMeshes.normal(::PhysicalMeshes.Polygon2D)
PhysicalMeshes.normal(::PhysicalMeshes.Polygon3D)
PhysicalMeshes.isconvex
PhysicalMeshes.circumcenter
```

### Rays

```@docs
PhysicalMeshes.intersect
PhysicalMeshes.reflect
```

## Mesh layer

### `MeshConfig` and traits

```@docs
PhysicalMeshes.MeshConfig
PhysicalMeshes.cell_volume
PhysicalMeshes.cell_area
PhysicalMeshes.cell_center
PhysicalMeshes.cell_neighbors
PhysicalMeshes.node_neighbors
PhysicalMeshes.mesh_dimension
PhysicalMeshes.number_of_nodes
PhysicalMeshes.number_of_cells
PhysicalMeshes.number_of_fields
PhysicalMeshes.add_node!
PhysicalMeshes.add_cell!
PhysicalMeshes.add_field!
PhysicalMeshes.get_field
PhysicalMeshes.interpolate_field
```

### Concrete meshes

```@docs
PhysicalMeshes.MeshCartesianStatic
PhysicalMeshes.DistributedMesh
PhysicalMeshes.calculate_local_config
PhysicalMeshes.partition_neighbors
```

### Field wrappers

```@docs
PhysicalMeshes.AbstractField
PhysicalMeshes.ScalarField
PhysicalMeshes.VectorField
PhysicalMeshes.TensorField
PhysicalMeshes.ArrayScalarField
PhysicalMeshes.ArrayVectorField
PhysicalMeshes.ArrayTensorField
```

### Mass-assignment operations

```@docs
PhysicalMeshes.assignmesh
PhysicalMeshes.mesh2particle
PhysicalMeshes.assignparticle
PhysicalMeshes.is_inbound
PhysicalMeshes.outbound_list
PhysicalMeshes.total_mass
PhysicalMeshes.mass_center
PhysicalMeshes.zoom
PhysicalMeshes.shrink
```

## Parallel layer

```@docs
PhysicalMeshes.ParallelBackend
PhysicalMeshes.detect_resources
PhysicalMeshes.select_backend
PhysicalMeshes.is_serial
PhysicalMeshes.is_threads
PhysicalMeshes.is_distributed
PhysicalMeshes.is_gpu
PhysicalMeshes.is_parallel
PhysicalMeshes.set_backend_fft_threads!
PhysicalMeshes.to_device
PhysicalMeshes.to_host
PhysicalMeshes.distribute
PhysicalMeshes.release!
PhysicalMeshes.parallel_sum
PhysicalMeshes.parallel_maximum
PhysicalMeshes.parallel_minimum
PhysicalMeshes.parallel_findmax
PhysicalMeshes.parallel_quantile
PhysicalMeshes.parallel_sumprod
PhysicalMeshes.parallel_broadcast!

PhysicalMeshes.exchange_boundary_data
PhysicalMeshes.sync_mesh
PhysicalMeshes.sync_field
```

## Algorithms

```@docs
PhysicalMeshes.AbstractFDTD
PhysicalMeshes.FDTD2D
PhysicalMeshes.FDTD3D
PhysicalMeshes.fdtd_step
PhysicalMeshes.fdtd_simulate
PhysicalMeshes.apply_boundary_conditions
PhysicalMeshes.apply_periodic_boundary
PhysicalMeshes.apply_dirichlet_boundary
PhysicalMeshes.apply_vacuum_boundary

PhysicalMeshes.AbstractPIC
PhysicalMeshes.PIC2D
PhysicalMeshes.PIC3D
PhysicalMeshes.pic_step
PhysicalMeshes.pic_simulate
PhysicalMeshes.boris_algorithm
PhysicalMeshes.interpolate_field
PhysicalMeshes.deposit_charge
PhysicalMeshes.deposit_current
PhysicalMeshes.relax_poisson
```

## Index

```@index
Pages = ["Methods.md"]
```

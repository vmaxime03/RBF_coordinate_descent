Smoother parameters file *Smoother_config.json*

Put the quad mesh and polyline in *assets/*

```json
{
  "bord":"bord.obj",                    // border for the smoother
  "quad":"quad_mesh.obj",               // quad mesh to smooth
  "smoother_w": double,                 // weight of the "border points stay on the border" constraint, need to be increased for lower target_interpolant_neighbors

  "nsample": int,                       // number of samples for each polyline segment 
  "target_function_width": double,      // distance field width around the polyline, depends of the mesh coordiante space
  "target_interpolant_neighbors": int,  // inter function influence (>= 2)
  "lambda_distance": double,            // weight of the 0 iso accuracy in the least square (1 is fine)
  "export_field": true                  // for distance fiel visualization
}
```


Build : 
```shell 
make build_release
```

Run smoother : 
```shell 
make run_smoother
```

Smoother output in *run_output/smoother/smoothed.obj*

Visualize ditance field (require Octave) :
```shell 
make colormap
```


### TODO 
- Refactor
- Find some of the parameters automaticaly (*smoother_w*, *taget_function_width*)
- Quickly find query point fonctions to sum over
- Multi frontier quad mmesh 



Smoother parameters file *Smoother_config.json*

Put the quad mesh and polyline in *assets/*

```json
{
  "bord":"bord.obj",                  // border for the smoother
  "quad":"quad_mesh.obj",             // quad mesh to smooth
  "smoother_w": 3,                    // weigth of the "border points stay on the border" constraint, work better when >= 3

  "nsample": 5,                       // number of samples for each polyline segment
  "target_function_width": 0.075,     // distance field width around the polyline
  "target_interpolant_neighbors": 6,  // inter function influence 
  "lambda_distance": 1,               // weigth of the 0 iso accuracy in the least square
  "export_field": true                // for distance fiel visualization
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


** TODO 
- Refactor
- Quickly find query point fonctions to sum over
- Multi frontier quad mmesh 



# Teeth Annotator

#### Introduction
* This code could be applied on normal size human gingival obj file.
* This code could differentiate each tooth from background gingival and each other. The cusp detection module could give predicted cusp location.

#### API
```
[F, X, frag_color] = annotate_teeth(filepath, anchor_points);
```
* Input: `anchor_points` is a three-elements array which should containe 3 points used for to locating a middle plane of the gingival.
* Output: `F` are fragments (each fragment consists of three IDs of verticles ), `X` are IDs and coordinates of verticles. `frag_color` is the classification array of every fragments, 0 for background and 1-n for n teeth.

#### Usage
* Before using, make sure to do:
```
mex build_graph.cpp utils.cpp;
mex buildE.cpp utils.cpp;
mex dijsktra.cpp utils.cpp
mex dye.cpp utils.cpp;
mex selectV.cpp utils.cpp;
```

#### Demo

![Teeth Annotation](/demo1.png)
![Cusp Detection](/demo2.png)

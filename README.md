simplicial-layout
=================
Code for drawing simplicial complexes in arbitrary dimensions.  Basically a generalized version of graph drawing for meshes and tetrahedra which works in 2D, 3D, and higher dimensions.

Usage
=====
First install the library like this:

    npm install simplicial-layout
    
Then you can create a graph layout as follows:

```javascript
var MeshLayout = require("simplicial-layout")

var l = new MeshLayout([[0, 1, 2], [2, 3], [3,4,5]], 2)
l.solve()
console.log(l.positions)
```

API
===

### `var MeshLayout = require("simplicial-layout")(cells[, dimension, options])`

* `cells` - A simplicial complex
* `dimension` - The dimension of the space to embed the complex in (default 2)
* `options` - Extra options for the layout.  Currently supports the following:
  + `radius` - Radius for each edge.  Default 2
  
### `layout.step(dt)`
Steps the mesh solver a finite amount

### `layout.solve()`
Solves the mesh layout completely

### `layout.positions`
Positions of the vertices

Credits
=======
(c) 2013 Mikola Lysenko. BSD



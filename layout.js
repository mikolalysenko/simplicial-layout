"use strict"

var EPSILON = 1e-4

var sc = require("simplicial-complex")
var numeric = require("numeric")

function compareLex(a,b) {
  for(var i=0, len=a.length-1; i<len; ++i) {
    var d = a[i] - b[i]
    if(d) { return d }
  }
  return 0
}

function Layout(cells, dimension, options) {
  dimension = dimension || 2
  options = options || {}

  var num_vertices = sc.countVertices(cells)
  
  this.cells      = sc.unique(sc.normalize(sc.explode(cells)))
  this.dimension = dimension|0
  this.numVertices = num_vertices|0
  this.radius     = options.radius || 1.0
  
  this.positions = numeric.rep([num_vertices, dimension], 0.0)
  this.forces = numeric.rep([num_vertices, dimension], 0.0)
  this.grid = numeric.rep([num_vertices<<dimension, dimension+1], 0)
  
  for(var i=0; i<num_vertices; ++i) {
    for(var j=0; j<dimension; ++j) {
      this.positions[i][j] = (Math.random() - 0.5) * num_vertices * this.radius
    }
  }
}

Layout.prototype.step = function(dt) {
  var dimension = this.dimension
  var radius = this.radius
  var logr = Math.log(this.radius)

  //Apply per-cell forces
  var centroid = new Array(dimension)
  var dir = new Array(dimension)
  for(var i=0, ilen=this.cells.length; i<ilen; ++i) {
    var c = this.cells[i]
    
    if(c.length <= 1) {
      continue
    }
    
    //Compute centroid
    for(var j=0; j<dimension; ++j) {
      centroid[j] = 0.0
    }
    for(var k=0, klen=c.length; k<klen; ++k) {
      var p = this.positions[c[k]]
      for(var j=0; j<dimension; ++j) {
        centroid[j] += p[j]
      }
    }
    var w = 1.0 / c.length
    for(var j=0; j<dimension; ++j) {
      centroid[j] *= w
    }
    
    //Apply forces
    for(var k=0, klen=c.length; k<klen; ++k) {
      var p = this.positions[c[k]]
      var f = this.forces[c[k]]
      var d = 0.0
      for(var j=0; j<dimension; ++j) {
        dir[j] = centroid[j] - p[j]
        d += dir[j] * dir[j]
      }
      var ds = Math.sqrt(d)
      if(Math.abs(ds) > EPSILON) {
        var m = (Math.log(ds) - logr) / ds
        for(var j=0; j<dimension; ++j) {
          f[j] += m * dir[j]
        }
      } else {
        for(var j=0; j<dimension; ++j) {
          f[j] += Math.random() - 0.5
        }
      }
    }
  }
  
  //Compute overlaps
  var coord = new Array(dimension)
    , ptr = 0
  for(var i=0, ilen=this.numVertices; i<ilen; ++i) {
    var p = this.positions[i]
    for(var j=0; j<dimension; ++j) {
      coord[j] = Math.floor(p[j] / radius)|0
    }
    for(var s=0, sl=1<<dimension; s<sl; ++s) {
      var g = this.grid[ptr++]
      for(var j=0; j<dimension; ++j) {
        g[j] = coord[j] + ((s>>>j)&1)
      }
      g[dimension] = i
    }
  }
  this.grid.sort(compareLex)

  //Process overlapping cells
  for(var i=0, len=this.grid.length; i<len;) {
    var cs = this.grid[i]
    var j=i
    while(++j < len && compareLex(this.grid[j], cs) === 0) {
      var a = this.grid[j][dimension]
      for(var k=i; k<j; ++k) {
        var b = this.grid[k][dimension]
        
        //FIXME: Check if this cell is the lexicographically smallest, if not skip
        
        //Apply interaction force to a/b
        var pa = this.positions[a]
        var fa = this.forces[a]
        var pb = this.positions[b]
        var fb = this.forces[b]
        
        var d = 0.0
        for(var l=0; l<dimension; ++l) {
          dir[l] = pa[l] - pb[l]
          d += dir[l] * dir[l]
        }
        var ds = Math.sqrt(d)
        if(ds < EPSILON) {
          for(var l=0; l<dimension; ++l) {
            fa[l] += Math.random() - 0.5
            fb[l] += Math.random() - 0.5
          }
        } else {
          var m = 10.0 / (ds * (1.0 + Math.exp(ds - 0.5*radius)))
          for(var l=0; l<dimension; ++l) {
            var f = m * dir[l]
            fa[l] += f
            fb[l] -= f
          }
        }
      }
    }
    i = j
  }
  
  //Move
  var tf = 0.0
  for(var i=0, ilen=this.positions.length; i<ilen; ++i) {
    var p = this.positions[i]
    var f = this.forces[i]
    for(var j=0; j<dimension; ++j) {
      p[j] += f[j] * dt
      tf = Math.max(tf, Math.abs(f[j]))
      f[j] = 0.0
    }
  }
  return tf
}

Layout.prototype.solve = function() {
  for(var i=0; i<this.numVertices*100; ++i) {
    if(this.step(0.1) > EPSILON) {
      break
    }
  }
}

module.exports = Layout
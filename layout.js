"use strict"

var EPSILON = 1e-4

var sc = require("simplicial-complex")
var numeric = require("numeric")
var uniq = require("uniq")

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

  this.lo = options.lo || numeric.rep([dimension], -Infinity)
  this.hi = options.hi || numeric.rep([dimension], Infinity)

  var num_vertices = sc.countVertices(cells)
  
  this.cells      = sc.unique(sc.normalize(sc.explode(cells)))
  this.stars      = sc.dual(cells, num_vertices)
  this.temperature = 1.0
  
  for(var i=0; i<this.stars.length; ++i) {
    var s = this.stars[i]
    var n = []
    for(var j=0; j<s.length; ++j) {
      n = n.concat(cells[s[j]])
    }
    n.sort()
    this.stars[i] = uniq(n)
  }
  
  this.dimension = dimension|0
  this.numVertices = num_vertices|0
  this.radius     = options.radius || 1.0
  
  this.positions = numeric.rep([num_vertices, dimension], 0.0)
  this.velocities = numeric.rep([num_vertices, dimension], 0.0)
  this.forces = numeric.rep([num_vertices, dimension], 0.0)
  this.grid = numeric.rep([num_vertices<<dimension, dimension+1], 0)
  
  for(var i=0; i<num_vertices; ++i) {
    if(options.lo && options.hi) {
      for(var j=0; j<dimension; ++j) {
        this.positions[i][j] = (Math.random() * (this.hi[j]-this.lo[j])) + this.lo[j]
      }
    } else {
      for(var j=0; j<dimension; ++j) {
        this.positions[i][j] = (Math.random() - 0.5) * num_vertices * this.radius
      }
    }
  }
}

Layout.prototype.jiggle = function() {
  this.temperature = 1.0
  for(var i=0; i<this.positions.length; ++i) {
    for(var j=0; j<this.dimension; ++j) {
      this.positions[i][j] += this.radius * (Math.random() - 0.5)
    }
  }
}

Layout.prototype.step = function(dt) {
  if(!dt) {
    dt = 0.1
  }
  var dimension = this.dimension
  var radius = this.radius
  var logr = Math.log(this.radius)

  //Apply per-cell forces
  var centroid = new Array(dimension)
  var dir = new Array(dimension)
  var cattract = radius
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
        var m = cattract * (Math.log(ds) - logr) / ds
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
    , repel_radius = 2.0*radius
  for(var i=0, ilen=this.numVertices; i<ilen; ++i) {
    var p = this.positions[i]
    for(var j=0; j<dimension; ++j) {
      coord[j] = Math.floor(p[j] / repel_radius)|0
    }
    for(var s=0, sl=1<<dimension; s<sl; ++s) {
      var g = this.grid[ptr++]
      for(var j=0; j<dimension; ++j) {
        g[j] = (coord[j] + ((s>>>j)&1))|0
      }
      g[dimension] = i|0
    }
  }
  this.grid.sort(compareLex)

  //Process overlapping cells
  var crepulse = (this.temperature * Math.pow(radius,dimension) + Math.pow(radius, 1.0/dimension))
  for(var i=0, len=this.grid.length; i<len;) {
    var cs = this.grid[i]
    var j=i
    
    while(++j < len && compareLex(this.grid[j], cs) === 0) {
      var a = this.grid[j][dimension]
      var na = this.stars[a]
k_loop:
      for(var k=i; k<j; ++k) {
        var b = this.grid[k][dimension]
        for(var l=0, ln=na.length; l<ln; ++l) {
          if(na[l] === b) {
            continue k_loop
          }
        }
        
        //Apply interaction force to a/b
        var pa = this.positions[a]
        var pb = this.positions[b]
        
        //Check if the coordinate is lexicographically first intersection
        for(var l=0; l<dimension; ++l) {
          var ac = Math.floor(pa[l]/repel_radius)|0
          var bc = Math.floor(pb[l]/repel_radius)|0
          if(ac > bc) {
            if(cs[l] !== ac) {
              continue k_loop
            }
          } else if(cs[l] !== bc) {
            continue k_loop
          }
        }
        
        //If so, then apply force
        var fa = this.forces[a]
        var fb = this.forces[b]
        
        var d = 0.0
        for(var l=0; l<dimension; ++l) {
          dir[l] = pa[l] - pb[l]
          d += dir[l] * dir[l]
        }
        if(d < EPSILON) {
          for(var l=0; l<dimension; ++l) {
            fa[l] += Math.random() - 0.5
            fb[l] += Math.random() - 0.5
          }
        } else {
          var ds = Math.sqrt(d)
          var m = crepulse / (ds * (1.0 + Math.exp((4.0-3.0*temperature)*(0.5*repel_radius-ds))))
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
    , lo = this.lo
    , hi = this.hi
  for(var i=0, ilen=this.positions.length; i<ilen; ++i) {
    var p = this.positions[i]
    var v = this.velocities[i]
    var f = this.forces[i]
    for(var j=0; j<dimension; ++j) {
      p[j] += v[j] * dt
      v[j] += f[j] * dt
      v[j] *= 0.9
      f[j] = 0.0
      tf = Math.max(tf, Math.abs(f[j]))
      
      if(p[j] < lo[j]) {
        p[j] = lo[j]
      } else if(p[j] > hi[j]) {
        p[j] = hi[j]
      }
    }
  }
  this.temperature = Math.max(0.0, this.temperature-dt/this.numVertices)
  
  return tf
}

Layout.prototype.solve = function() {
  for(var i=0; i<this.numVertices*100; ++i) {
    var f = this.step(0.1)
    if(f < EPSILON) {
      break
    }
  }
}

module.exports = Layout
var MeshLayout = require("../layout.js")

require("tap").test("checking layout", function(t) {

  //var cells = [[0,1,2],[3],[4], [5,6]]
  //var cells = [[0,1,2,3]]
  var cells = [[0,1,2],[2,3,6],[1,4,7],[3,4,5]]
  
  
  var l = new MeshLayout(cells, 2)
  
  l.solve()
  console.log(l.positions)

  t.end()
})
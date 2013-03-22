var createLayout = require("../layout.js")

require("tap").test("checking layout", function(t) {

  var cells = [[0,1], [1, 2], [2, 3]]
  
  var l = createLayout(cells, 2)
  
  l.solve()
  console.log(l.positions)
  


  t.end()
})
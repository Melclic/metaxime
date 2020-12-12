var radius = 100;

var fill = d3.scale.category20();

var force = d3.layout.force()
    .gravity(.1)
    .charge(-1800)
    .linkDistance(350)
    .linkStrength(0.001)
    .size([document.body.clientWidth, document.body.clientHeight])
    .friction(0.9)
    .theta(.8)
    .alpha(.1)

var svg = d3.select("body")
    .append("svg")
      .attr("width", "100%")
      .attr("height", "100%")
      .call(d3.behavior.zoom().on("zoom", function () {
	  svg.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
      }))
    .append("g");


d3.

d3.json("graph.json", function(error, graph) {
  if (error) throw error;

  var link = svg.selectAll("line")
      .data(graph.links)
      .enter().append("line");

  /*
  var node = svg.selectAll("circle")
      .data(graph.nodes)
      .enter().append("circle")
      .attr("r", radius - .75)
      .style("fill", function(d) { return fill(d.group); })
      .style("stroke", function(d) { return d3.rgb(fill(d.group)).darker(); })
      .call(force.drag);
  */

  var node = svg.selectAll("rect")
      .data(graph.nodes)
      .enter().append("rect")
      .attr("width", radius)
      .attr("height", radius)
      .style("fill", function(d) { return fill(d.group); })
      .style("stroke", function(d) { return d3.rgb(fill(d.group)).darker(); })
      .call(force.drag);



  makeStrc('[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]', 'example-strc');
  const svg1 = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg1.setAttribute("width", "100");
  svg1.setAttribute("height", "100");
  //makeStrc('[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]', svg1);
  console.log(svg1)

  force.nodes(graph.nodes)
       .links(graph.links)
       .on("tick", tick)
       .start();

  function makeStrc(input, canvas) {
    smilesDrawer = new SmilesDrawer.Drawer({});
    SmilesDrawer.parse(input, function(tree) {
          // Draw to the canvas
          smilesDrawer.draw(tree, canvas, "light", false);
          // Alternatively, draw to SVG:
          //svgDrawer.draw(tree, canvas, 'dark', false);
        });
  }

  console.log(svg.select('example-strc'))

  function tick() {

    node.attr("x", function(d) {return d.x = Math.max(radius, Math.min( document.body.clientWidth - radius, d.x)); })
        .attr("y", function(d) {return d.y = Math.max(radius, Math.min( document.body.clientHeight - radius, d.y)); });

    node.each(function(d){w=200*(1+d.group); d.x -= (.2* (d.x-w)) })

    link.attr("x1", function(d) { return d.source.x+radius/2; })
        .attr("y1", function(d) { return d.source.y+radius/2; })
        .attr("x2", function(d) { return d.target.x+radius/2; })
        .attr("y2", function(d) { return d.target.y+radius/2; });

  }
});

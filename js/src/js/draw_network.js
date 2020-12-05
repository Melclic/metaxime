function drawNetwork(svg, network_file, show_only_central=true, chem_width=100, chem_height=100) { 
  // Create the input graph
  var g = new dagreD3.graphlib.Graph()
    .setGraph({nodesep: 30,
	       ranksep: 150,
	       rankdir: "LR",
	       marginx: 10,
	       marginy: 10})
    .setDefaultEdgeLabel(function() { return {}; });
  d3.json(network_file).then(function(graph) {
    console.log(graph)
    var ignore_species = [];
    //################### nodes ##################
    for (var i = 0; i < graph.nodes.length; i++) {
      if (show_only_central) { 
	if (!graph.nodes[i].central_species && graph.nodes[i].type=='species') {
	  ignore_species.push(graph.nodes[i].id);
	  continue;
	}
      }
      if (graph.nodes[i].type=='species') {
	var chem_svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
	chem_svg.id = 'svg_'+graph.nodes[i].id;
	chem_svg.setAttribute("width", chem_width);
	chem_svg.setAttribute("height", chem_height);
	//var serializer = new XMLSerializer();
	//var source = serializer.serializeToString(chem_svg);
	g.setNode(graph.nodes[i].id, {labelType: 'html', //can also use 'svg' but not centered
				      label: chem_svg,
				      width: chem_width,
				      height: chem_height}
	);
      }
      else if(graph.nodes[i].type=='reaction') {
	g.setNode(graph.nodes[i].id, {labelType: 'string', 
				      label: graph.nodes[i].id,
				      height: chem_height/2,
				      width: chem_width}); 
      }
      else {
	console.log('ERROR: Cannot interpret the following node type: '+graph.nodes[i].type);
      }
    }
    //################### links ##################
    for (var i = 0; i < graph.links.length; i++) {
      if(!ignore_species.includes(graph.links[i].source) && !ignore_species.includes(graph.links[i].target)) {
	g.setEdge(graph.links[i].source, graph.links[i].target);
      }
    }
    g.nodes().forEach(function(v) {
      var node = g.node(v);
      // Round the corners of the nodes
      node.rx = node.ry = 5;
    });
    // Create the renderer
    var render = new dagreD3.render();
    // Set up an SVG group so that we can translate the final graph.
	//WARNING: This overwrites the passed svg HTMLElement
    //var svg = d3.select("svg");
    var svgGroup = svg.append("g");
    // Run the renderer. This is what draws the final graph.
    render(d3.select("svg g"), g);
    // Center the graph
    //var inner = svg.firstChild;
    //var xCenterOffset = (svg.attr("width") - g.graph().width) / 2;
    //inner.attr("transform", "translate(" + xCenterOffset + ", 20)");
    svg.attr("height", g.graph().height + 40);

    
    //svg.attr("height", g.graph().height + 40);
    //draw the molecules in the svg space
    for (var i = 0; i < graph.nodes.length; i++) {
      if (graph.nodes[i].type=='species' && !ignore_species.includes(graph.nodes[i].id)) {
	canvas_static = document.getElementById('svg_'+graph.nodes[i].id)
	var explicitH = false;
	if (graph.nodes[i].id=='MNXM1__64__MNXC3') {
	  explicitH = true;
	}
    console.log(graph.nodes[i].brsynth.smiles);
	//makeSVGstrc(graph.nodes[i].brsynth.smiles.toUpperCase(), //WARNING: lower case means aromatic atoms. 
	makeSVGstrc(graph.nodes[i].brsynth.smiles, //WARNING: lower case means aromatic atoms. 
		    'svg_'+graph.nodes[i].id,
		    chem_width,
		    chem_height,
		    explicitH) 
      }
    }
  });
  
  d3.selectAll('node')
    .attr('test', function(d) { console.log(d) })
    .call(d3.drag().on('drag', dragged));


  //################ functions ########################
  function dragged(d) {
    d.x = d3.event.x, d.y = d3.event.y;
    d3.select(this).attr("cx", d.x).attr("cy", d.y);
    link.filter(function(l) { return l.source === d; }).attr("x1", d.x).attr("y1", d.y);
    link.filter(function(l) { return l.target === d; }).attr("x2", d.x).attr("y2", d.y);
  }
  function makeSVGstrc(input, in_html_element, in_width, in_height, explicitH) {
    svgDrawer = new SmilesDrawer.SvgDrawer({width: in_width,
					    height: in_height,
					    explicitHydrogens: explicitH});
    SmilesDrawer.parse(input, function(tree) {
      svgDrawer.draw(tree, in_html_element, 'light', false);
    });
  }
}

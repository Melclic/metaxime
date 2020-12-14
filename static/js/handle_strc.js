//Function that uses the smilesDrawer to draw in canvas a smiles structure
function makeStrc(input, canvas, width, height, bck_color='light') {
  var smilesDrawer = new SmilesDrawer.Drawer({ width: width, height: height });
  SmilesDrawer.parse(input, function(tree) {
    smilesDrawer.draw(tree, canvas, bck_color, false);
  });
  //TODO: need to send error message when invalid input
}

//Function that quesries a web service to retrun
function searchSMILESbyName(in_name) {
  var smiles = null;
  if (in_name) {
    $.ajax({
      url: 'https://cactus.nci.nih.gov/chemical/structure/'+in_name+'/smiles',
      async: false,
      success: function(result) {
	if (result) {
	  smiles = result;
	}
	else {
	  console.log('searchSMILESbyName: Empty string returned');
	  return '';
	}
      },
      error: function(xhr) {
	//alert("An error occured: " + xhr.status + " " + xhr.statusText);
	console.log('searchSMILESbyName: Error message from API call');
	console.log(xhr.status);
	return '';
      }
    });
  }
  else {
    console.log('searchSMILESbyName: Invalid entry');
    return '';
  }
  return smiles;
}


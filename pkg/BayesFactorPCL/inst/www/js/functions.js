function setup(){
	$('#logTH').prepend( textBase );
}

function listBayesFactors(bayesFactors){
	var i;
	var n = bayesFactors.length;
	var checked;
	$('#bfTableBody').html('');
	
	for(i = 0; i < n; i++){
		if(bayesFactors[i].isBase){
			checked = "checked='checked'";
		}else{
			checked = "";
		}
		checkBox = "<input name='baseM' type='radio' " + 
						checked + 
						" value='" + bayesFactors[i].model + "' " + 
						"/>";
		$('#bfTableBody').append('<tr>' +
										  '<td>' + bayesFactors[i].model + '</td>' +
										  '<td>' + checkBox + '</td>' +
										  '<td>' + bayesFactors[i].name + '</td>' +
										  '<td>' + bayesFactors[i].bf + '</td>' +
										  '<td>' + Math.log(bayesFactors[i].bf) / Math.log(base) + '</td>' +
										  '<td>' + bayesFactors[i].iterations + '</td>' +
										  '</tr>'
										  );
	}
	$("input[name=baseM]").change( changeBase );
}

function changeBase(){
	//alert(this.value);
}

function listEffects(){
	$('#effectsTableBody').html('');
	var checkBox;
	var way;
	var nFac;
	var i;
	var wayText;
	
	$.getJSON("/custom/aov/data?what=nFac", 
		function(data) { 
			nFac = parseInt(data);
			for(i=0;i<nFac;i++){
				if(i==0){
					wayText = "Main effects";
				}else{
					wayText = (i+1) + " way";
				}
				$('#effectsTableBody').append(
 					"<tr id='" + (i+1) + "way'>" +
 					'<td>' + wayText + '</td>' + 
 					'<td></td>' +
 					'<td></td>' +
 					'</tr>'
 				);	 			
	
				$("#" + (i+1) + "way").click( function(){
					$(".row" + this.id).toggle("fast");
				});
			}
		});
	 
	$.getJSON("/custom/aov/data?what=fixed",
  		function(data) {
 			$.each(data, function(index,value){
				checkBox = "<input id='effect" + (index+1) + "' type='checkbox' " + "/>";
 				way = value.split(":").length;
 				$('#' + way + "way").after(
 					"<tr class='row" + way + "way'>" +
 					'<td>' + (index+1) + '</td>' + 
 					'<td>' + checkBox + '</td>' +
 					'<td>' + value + '</td>' +
 					'</tr>'
 				);	 			
 			});
 		});

	for(i=0;i<nFac;i++){

	}

}

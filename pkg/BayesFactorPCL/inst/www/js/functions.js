function setup(){
	$('#logTH').prepend( textBase );
	
	$('#analyzeSelectedButton').click( function(){
		var whichModel = chosenModel();
		analyzeModel( whichModel ); 
	});
}

function listBayesFactors(){
	var i;
	var n = bayesFactors.length;
	var checked;
	$('#bfTableBody').html('');
	var baseBF = baseBayesFactor();
	
	for(i = 0; i < n; i++){
		if(bayesFactors[i].isBase){
			checked = "checked='checked'";
		}else{
			checked = "";
		}
		checkBox = "<input name='baseM' type='radio' " + 
						checked + 
						" value='" + i + "' " + 
						"/>";
		$('#bfTableBody').append('<tr>' +
										  '<td>' + bayesFactors[i].model + '</td>' +
										  '<td>' + checkBox + '</td>' +
										  '<td>' + bayesFactors[i].name + '</td>' +
										  '<td>' + Math.exp( bayesFactors[i].bf - baseBF ) + '</td>' +
										  '<td>' + (bayesFactors[i].bf  / Math.log(base) - 
										  		    baseBF  / Math.log(base)) + '</td>' +
										  '<td>' + bayesFactors[i].iterations + '</td>' +
										  '<td>' + bayesFactors[i].rscaleFixed + '</td>' +
										  '<td>' + bayesFactors[i].rscaleRandom + '</td>' +
										  '</tr>'
										  );
	}
	$("input[name=baseM]").change( changeBase );
}

function changeBase(){
	var i;
	for(i=0;i<bayesFactors.length;i++){
		if(this.value == i) {
			bayesFactors[i].isBase = true;
		}else{
			bayesFactors[i].isBase = false;
		}
	}
	makePlot();
	listBayesFactors();
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
				checkBox = "<input id='effect:" + index + "' type='checkbox' " + "/>";
 				way = value.split(":").length;
 				$('#' + way + "way").after(
 					"<tr class='row" + way + "way'>" +
 					'<td>' + index + '</td>' + 
 					'<td>' + checkBox + '</td>' +
 					'<td>' + value + '</td>' +
 					'</tr>'
 				);	 			
 			});
 		});

}

function chosenModel(){
	var whichModel = 0;
	var effNum;
	$('input:checked[id^="effect:"]').each( function() {
		effNum = parseInt(this.id.split(':')[1]);
		whichModel += Math.pow(2, effNum);
	});
	return(whichModel);
}

function analyzeModel(whichModel) {
	var iterations = $("#nIterations").val();
	var rscaleFixed = $("#rscaleFixed").val();
	var rscaleRandom = $("#rscaleRandom").val();
	
	$.getJSON("/custom/aov/data?", 
		{
			what: "analysis",
			model: whichModel,
			iterations: iterations,
			rscaleFixed: rscaleFixed,
			rscaleRandom: rscaleRandom
		},
		setResults );
}

function setResults(data) {
	bayesFactors.push(data);
	if(bayesFactors.length == 1){
		bayesFactors[0].isBase = true;
	}else if(bayesFactors.length > 1){
		makePlot();
	}
	listBayesFactors();
}

function baseBayesFactor(){
	var whichBase; 
	var i;
	for(i=0;i<bayesFactors.length;i++){
		if(bayesFactors[i].isBase){
			whichBase = i;
		}
	}
	return(bayesFactors[whichBase].bf);
}

function makePlot()
{
	var baseBF = baseBayesFactor();
	$("#bfImageContainer").html("<img src='/custom/aov/bf.png?baseBF=" + baseBF + "'/>");
}
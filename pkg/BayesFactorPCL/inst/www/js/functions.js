
// http://stackoverflow.com/questions/6967975/how-to-generate-and-append-a-random-string-using-jquery
function randString(n)
{
    if(!n)
    {
        n = 5;
    }

    var text = '';
    var possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

    for(var i=0; i < n; i++)
    {
        text += possible.charAt(Math.floor(Math.random() * possible.length));
    }

    return text;
}

// http://stackoverflow.com/questions/2595835/how-to-hide-table-columns-in-jquery
function hideExtraColumns() {
	$.each(advCols, function(index, value){
			$('#bfTable td:nth-child('+(value+1)+')').hide();
			$('#bfTable th:nth-child('+(value+1)+')').hide();
		});
}

function showExtraColumns() {
	$.each(advCols, function(index, value){
			$('#bfTable td:nth-child('+(value+1)+')').show();
			$('#bfTable th:nth-child('+(value+1)+')').show();
		});
}

function update(plot){
		plot = (typeof plot === "undefined") ? true : plot;
		if(plot) makePlot();
		makeSaveLink();
		listBayesFactors();
		changeProgress();
}

function changeLogBase(){
		textBase = $('#logBase option:selected').text();
		$('#logTHspan').html( textBase );
}

function setup(){
	setDefaults();
	changeLogBase();
	$("#logBase").change( function() { changeLogBase(); update(); } );
	
	$('#analyzeSelectedButton').click( function(){
		var whichModel = chosenModel();
		analyzeModels( [ whichModel ] ); 
	});
	$('#analyzeAllButton').click( allNways );
	$('#analyzeTopButton').click( topNways );
	$('#clearButton').click( clearAnalyses );
	$('#cullButton').click( cullAnalyses );

	
	$('#bfTableContainer').jScrollPane();
	$('#effectsTableContainer').jScrollPane();
	$('#extraColsToggle').change( checkExtraCols );
	
	$('#bfImageContainer').click( makePlot );
	
	$("#showdiv").click( toggleOptions );
	
	$(".effectsSelect").click( effectsSelect );
	
	roundTo = parseInt($("#roundTo").val());
	$("#roundTo").change( function() { 
		roundTo = parseInt($(this).val());
		listBayesFactors();
	});

	$.ajaxSetup({ cache:false });
}

function cullAnalyses(plot){
	var cullLevel = parseFloat($('#cullLevel').val());
	var culled = [];
	var baseBF = baseBayesFactor();
	jQuery.each(bayesFactors, function(index, value){
		var bf = value.bf - baseBF;
		if(bf > -Math.log(cullLevel)){
			culled.push(value);
		}
	});
	bayesFactors = culled;
	update(plot);
}

function autoCull(){
	if($("#autoCull").prop("checked")){
		cullAnalyses(false);
	}
}

function baseIsMax(){
	var max = Number.NEGATIVE_INFINITY;
	var whichMax = -1;
	jQuery.each(bayesFactors, function(index, value){
		if(value.bf > max){
			max = value.bf;
			whichMax = index;
		}
	});
	$('[name=baseM][value="' + whichMax +'"]').attr('checked', 'checked');
	$('[name=baseM][value="' + whichMax +'"]').trigger('change');
}

function setBase(){
	if($("#baseIsMax").prop("checked")){
		baseIsMax();
	}
}

function toggleOptions(){
	if ( $("#optionsContainer").is(':visible') ){
		$("#optionsContainer").slideUp(300);
		$("#optionsContainerIcon").html("[+]");
	}else{
		$("#optionsContainer").slideDown(300);
		$("#optionsContainerIcon").html("[-]");
	}

}

function effectsSelect(){
	if( this.id == "effectsSelectNone" ){
		$('#effectsTable input[id^="effect"]').prop('checked', false);
	}
	if( this.id == "effectsSelectAll" ){
		$('#effectsTable input[id^="effect"]').prop('checked', true);
	}	
}

function clearAnalyses(){
	bayesFactors = [];
	update();
}

function checkExtraCols(){
	var checked = $('#extraColsToggle').is(':checked');
	if(checked){
		showExtraColumns();
	}else{
		hideExtraColumns();
	}
}

function listBayesFactors( sortby ){
	var i;
	var n = bayesFactors.length;
	var checked;
	var base = $('#logBase').val();

	$('#bfTableBody').html('');
	if(n<1) {return;}
	var baseBF = baseBayesFactor();
	
	if(typeof sortby !== "undefined"){
		bayesFactors.sort(dynSort(sortby));
	}
	
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
		
		var bf, logbf;
		if(bayesFactors[i].bf == null){
			bf = "-";
			logbf = "-";			
		}else{
			bf = Math.exp( bayesFactors[i].bf - baseBF );
			logbf = (bayesFactors[i].bf - baseBF  ) / Math.log(base);
		}
		
		$('#bfTableBody').append('<tr id="bfRow-' + i +  '">' +
										  '<td>' + bayesFactors[i].model + '</td>' +
										  '<td>' + checkBox + '</td>' +
										  '<td class="status"><div class="progressBarContainer">' +
										  		'<div class="progressBar"></div></div>' + bayesFactors[i].status + '</td>' +
										  '<td id="niceName-' + i +  '">' + bayesFactors[i].niceName + '</td>' +
										  '<td id="name-' + i +  '" class="modelName">' + bayesFactors[i].name + '</td>' +
										  '<td class="roundMe">' + bf + '</td>' +
										  '<td class="roundMe">' + logbf + '</td>' +
										  '<td>' + bayesFactors[i].iterations + '</td>' +
										  '<td>' + bayesFactors[i].rscaleFixed + '</td>' +
										  '<td>' + bayesFactors[i].rscaleRandom + '</td>' +
										  '<td>' + bayesFactors[i].time + '</td>' +
										  '<td>' + bayesFactors[i].duration + '</td>' +
										  '<td id="refresh-' + i + '"><img src="img/refresh.png" class="refreshimg"/></td>' +
										  '<td id="delete-' + i + '">&#10006;</td>' +
										  '</tr>'
										  );

		if( i%2 ) {
			$("#bfTableBody > tr").last().addClass("odd");
		}
		if(bayesFactors[i].bf == null){
			$('[name=baseM][value="' + i +'"]').attr('disabled', 'disabled');
		}
				
	}
	$('td[id^="refresh-"]').click( refreshRow );
	$('td[id^="delete-"]').click( deleteRow );
	$('td[id^="delete-"]').addClass( "deleteButton" );		
	$('td[id^="name-"]').click( selectedModel );
	$('td[id^="niceName-"]').attr('contentEditable',true);
	$('td[id^="niceName-"]').blur( changeNiceName );
	$("input[name=baseM]").change( changeBase );
	$(".roundMe").each( roundCell );
	$('#bfTableContainer').data('jsp').reinitialise();
	checkExtraCols();
}

function selectedModel() {
	var i = this.id.split("-")[1];
	var mod = parseInt(bayesFactors[i].model);
	var bin = mod.toString(2);    
	var ar = bin.split("");
	var eff;
	$('input[id^="effect-"]').prop('checked', false);
	for(i=0;i<ar.length;i++){
    	eff = ar.length - 1 - i;
    	if(parseInt(ar[i])){
       	  $('#effect-' + eff).prop('checked', true);
    	}
    }
	
}




function changeProgress()
{
	$("td.status:contains('%')").each(function(index){
		var percent = parseInt($(this).text().split("%")[0]);
		$(this).children('.progressBarContainer').css("height","10px");
		$(this).children('.progressBarContainer').css("width","100%");
		$(this).children('.progressBarContainer').css("border","1px solid gray");
		$(this).find('.progressBar').html("&nbsp;");
		$(this).find('.progressBar').css("background-color","rgba(0,0,0,0.2)");
		$(this).find('.progressBar').css("height","100%");
		$(this).find('.progressBar').css("width",percent/100*40);
		//$(this).find('.progressBar').css("border","1px solid gray");
	});
}

function roundCell(){
	var value = parseFloat($(this).text());
	var pow10 = Math.pow(10,roundTo);
	$(this).html(Math.round(value*pow10)/pow10);
}

function changeNiceName() {
	var i = this.id.split("-")[1];
	bayesFactors[i].niceName = $(this).text();
	$(this).html($(this).text());
	makePlot();
}

function refreshRow() {
	var i = this.id.split("-")[1];
	var mod = bayesFactors[i].model;
	analyzeModels([ mod ], true, [ i ]);
}


function deleteRow() {
	var i = this.id.split("-")[1];
	var	baseBF = baseBayesFactor("index");
	bayesFactors.splice(i,1);
	if(bayesFactors.length>0 && baseBF==i){
			bayesFactors[0].isBase = true;
	}
	update();
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
	update();
}

function allNways(){
	$("#bfImageContainer").html("Click to refresh.");
	var i;
	if(nFac > 4) { return; }
	var toAnalyze = [];	
	
	var nModels = Math.pow(2, Math.pow(2, nFac) - 1);
	for(i=0;i<nModels;i++){
		toAnalyze.push(i);
	}
	analyzeModels( toAnalyze, true );
}

function topNways(){
	$("#bfImageContainer").html("Click to refresh.");
	
	var i, binaryRep;
	var nChar = Math.pow(2, nFac) - 1;
	var topModel = Math.pow(2, nChar) - 1;
	var toAnalyze = [];
	
	toAnalyze.push(topModel);
	
			
	for(i=0;i<nChar;i++){
		binaryRep = Array(nChar + 1).join("1").split("");
		binaryRep[i] = "0";
		binaryRep = binaryRep.join("");
		toAnalyze.push( parseInt(binaryRep,2));
	}		
	analyzeModels( toAnalyze, true );

}


function listEffects(){
	$('#effectsTableBody').html('');
	 
	$.getJSON("/custom/aov/data?what=fixed",
  		function(data) {
 			var effectNames = data.effects;
 			var checkBox;
			var way;
			var i;
			var wayText;
			
 			for(i=0;i<nFac;i++){
				if(i==0){
					wayText = "Main effects";
				}else{
					wayText = (i+1) + " way";
				}
				$('#effectsTableBody').append(
 					"<tr id='" + (i+1) + "way' class='effectsTableCategory'>" +
 					'<td><span class="effectsTableCatToggle">[-]</span> ' + wayText + '</td>' + 
 					'<td></td>' +
 					'<td></td>' +
 					'</tr>'
 				);	 			
			}
 					
 			$(".effectsTableCategory").click( function(){
				if($(this).find(".effectsTableCatToggle").html() == "[-]"){
					$(this).find(".effectsTableCatToggle").html("[+]")
				}else{
					$(this).find(".effectsTableCatToggle").html("[-]")
				}
				$(".effectsRow" + this.id).toggle("fast");
			});
			
 			$.each(effectNames, function(index,value){
				checkBox = "<input id='effect-" + index + "' type='checkbox' " + "/>";
 				way = value.split(":").length;
 				$('#' + way + "way").after(
 					"<tr class='effectsRow" + way + "way'>" +
 					'<td>&nbsp;</td>' + 
 					'<td>' + value + '</td>' +
 					'<td>' + checkBox + '</td>' +
 					'</tr>'
 				);	 			
 			});
 		});
	$('#effectsTableContainer').data('jsp').reinitialise();
}

function chosenModel(){
	var whichModel = 0;
	var effNum;
	$('input:checked[id^="effect-"]').each( function() {
		effNum = parseInt(this.id.split('-')[1]);
		whichModel += Math.pow(2, effNum);
	});
	return(whichModel);
}

function createTokens(n){
	var i;
	var tokens = [];
	for(i=0;i<n;i++){
		tokens.push("a" + randString(10));
	}
	return(tokens);
}

function analyzeModels(whichModels, plot, placeat) {
	plot = (typeof plot === "undefined") ? "true" : plot;
	placeat = (typeof placeat === "undefined") ? new Array(whichModels.length) : placeat;
	
	var tokens = createTokens(whichModels.length);
	var i;
	
	var iterations = $("#nIterations").val();
	var rscaleFixed = $("#rscaleFixed").val();
	var rscaleRandom = $("#rscaleRandom").val();
	
	for(i=0;i<whichModels.length;i++){
		addNewBayesFactor(tokens[i], whichModels[i], iterations, rscaleFixed, rscaleRandom, placeat[i]);
		update(false);
	}
	
	
	$.getJSON("/custom/aov/data?", 
		{
			what: "analysis",
			models: whichModels.join(","),
			iterations: iterations,
			rscaleFixed: rscaleFixed,
			rscaleRandom: rscaleRandom,
			tokens: tokens.join(",")
		},
		function(data) { startAnalysis(data, plot); });

	intervalRef = window.setInterval(function(){
		jQuery.getJSON("/custom/aov/update?", 
				{
					tokens: tokens.join(","),
					time: ((new Date())+"")
				}, function(data){ updateProgressHandler(data, plot); }
			);
	}, progressPollTime);
	
}

function updateProgressHandler(data, plot) {
	var tokens = data.tokens;
	var percents = data.percents;
	var statuses = data.statuses;
	var allstatus = data.allstatus;
	
	var status, percent, row, base;
		
	jQuery.each(tokens, function(index,token){
		status = statuses[index];
		percent = parseInt(percents[index]);
		row = tokenRow(token);		
		base = bayesFactors[row].isBase;
		var oldNiceName = bayesFactors[row].niceName;
		bayesFactors[row] = jQuery.parseJSON(data.returnLists[index]);
		
		// replace nice name with previous (for refreshes)
		if(oldNiceName !== undefined) bayesFactors[row].niceName = oldNiceName;
		
		bayesFactors[row].isBase = base;
		if(status=="done"){
			bayesFactors[tokenRow(token)].status = "";
		}else if(status=="running"){
			bayesFactors[tokenRow(token)].status = percent + "%";
		}
	});
	
	if(allstatus=="done"){
		window.clearInterval(intervalRef);
		setBase();
		autoCull();
		update(plot);
	}else{
		update(false);
	}
}

function addNewBayesFactor(token, whichModel, iterations, rscaleFixed, rscaleRandom, placeat){
	if( ( typeof placeat === "undefined" )  || ( placeat == null ) ){
		bayesFactors.push({
			model: whichModel,
			iterations: iterations,
			rscaleFixed: rscaleFixed,
			rscaleRandom: rscaleRandom,
			token: token,
			isBase: false,
			status: "0%"
		});
	}else if( placeat < bayesFactors.length ){
		bayesFactors[placeat] = 
		{
			model: whichModel,
			name: bayesFactors[placeat].name,
			niceName: bayesFactors[placeat].niceName,
			iterations: iterations,
			rscaleFixed: rscaleFixed,
			rscaleRandom: rscaleRandom,
			token: token,
			isBase:  bayesFactors[placeat].isBase,
			status: "0%"
		} 
	}
	if(bayesFactors.length == 1){
		bayesFactors[0].isBase = true;
	}
}

function startAnalysis(data, plot) {
	if(data.status=="busy"){
		alert("Busy response from BayesFactor.");
	}else if(data.status=="started"){
		//alert("Analyses started.");
	}else{
		alert("Invalid response from BayesFactor.");
	}
}


function baseBayesFactor(which){
	var whichBase=-1; 
	var i;
	for(i=0;i<bayesFactors.length;i++){
		if(bayesFactors[i].isBase){
			whichBase = i;
			break;
		}
	}
	if(which=="index"){
		return(whichBase);
	}else if (whichBase == -1){
		return(0);
	}else{
		return(bayesFactors[whichBase].bf);
	}
}

function tokenRow(token){
	var which; 
	var i;
	for(i=0;i<bayesFactors.length;i++){
		if(bayesFactors[i].token == token){
			return(i);
		}
	}
	return(null);
}


function extractBayesFactors(){
	allBFs = [];
	var i;
	var baseBF = baseBayesFactor();
	for(i=0;i<bayesFactors.length;i++){
			allBFs.push( 
				[ Math.exp( (bayesFactors[i].bf - baseBF) ), i ]
			);
	}
	return(allBFs);
}

function extractNames(){
	allNames = [];
	var i;
	var baseBF = baseBayesFactor();
	for(i=0;i<bayesFactors.length;i++){
			allNames.push( bayesFactors[i].name );
	}
	return(allNames);
}

function makePlot()
{
	var base = $('#logBase option:selected').text();
	var culled = [];
	jQuery.each(bayesFactors, function(index, value){
		if(value.bf != null){
			culled.push(value);
		}
	});
	if(culled.length>1){
		var BFobj = JSON.stringify(culled);
		var qstr = jQuery.param({logBase: base, BFobj: BFobj })
		$("#bfImageContainer").html("<img src='/custom/aov/bfs.png?" + qstr + "'/>");
	}else{
		$("#bfImageContainer").html("Analyze models for plot.");
	}
}

function makeSaveLink()
{
	var base = $('#logBase option:selected').text();
	var culled = [];
	jQuery.each(bayesFactors, function(index, value){
		if(value.bf != null){
			culled.push(value);
		}
	});
	if(culled.length>0){
		var BFobj = JSON.stringify(culled);
		var qstr = jQuery.param({logBase: base, BFobj: BFobj })
		$("#bfCSVlink").html("<a href='/custom/aov/saveobj?" + qstr + "&which=CSV'>CSV download</a>");
	}else{
		$("#bfCSVlink").html("");
	}
}


function dynSort(property) {
    return function (a,b) {
        return (a[property] < b[property]) ? -1 : (a[property] > b[property]) ? 1 : 0;
    }
}

function setDefaults(){
	jQuery.getJSON('/custom/aov/getdefaults', function(data){
		$("#rscaleFixed").val(data.rscaleFixed);
		$("#rscaleRandom").val(data.rscaleFixed);
		setDefaultIterations(data.iterations);
	});
}

function setDefaultIterations(iterations){
		var found = false;
		$("#nIterations").children("option").each(function(){
    		if(parseInt($(this).attr("value")) == iterations){
    			found = true;
    		}
		});
		if(!found){
			$("#nIterations").prepend(
				"<option value='" +iterations + "'>" + numberWithCommas(iterations) + "</option>"
			);
		}
		$("#nIterations").val(iterations);
}

// From http://stackoverflow.com/questions/2901102/how-to-print-number-with-commas-as-thousands-separators-in-javascript
function numberWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}
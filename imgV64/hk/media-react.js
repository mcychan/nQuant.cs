<!--
var root = 'media';
var cellsPerRow = 5;
var items;

function getItem(index)
{     
	var item = {};
    try
    {
		item.url = getElementTextNS("", "url", items[index], 0);
        item.title = getElementTextNS("", "title", items[index], 0);
        item.link = getElementTextNS("", "link", items[index], 0);
    }
    catch(err)
    {
        item.title = '&nbsp;';
        item.link = '#';
    }
	return item;
}

function renderHeader(header, tbId)
{	
	if(header.trim().length <=0)
		return "";
	
	var tbody = document.getElementById("tbMain");
	if(tbody && tbody.tagName == 'DIV')
		return "<div class='col'><h2>" + header + "</h2></div>";

	tbody = document.createElement('tbody');
	var tr = tbody.insertRow(tbody.rows.length);
    var td = tr.insertCell(tr.cells.length);
	td.innerHTML = "<b style='color: #00007D'>" + header + "</b>";

	tr = tbody.insertRow(tbody.rows.length);
	td = tr.insertCell(tr.cells.length);
	td.innerHTML = "<table id='" + tbId + "' bgColor='#808080' style='display:block' width='100%' cellSpacing='1' cellPadding='3'></table><hr />";
	return tbody.innerHTML;
}

function renderPage(doc)
{	
	var types = doc.getElementsByTagName(root)[0].childNodes;
	var resultHtml = "";

	for(var i=0; i<types.length; ++i)
	{	
		if(types[i].nodeName == '#text')
			continue;	
		resultHtml += renderHeader(types[i].getAttribute("text"), types[i].nodeName);
		items = doc.getElementsByTagName(types[i].nodeName)[0].childNodes;
		resultHtml += renderTables(types[i].nodeName);
	}	
	return resultHtml;
}

function renderTables(tbId)
{
	var tbody = document.getElementById(tbId);
	if(!tbody)
		return "";
	
	var tmp = 0;
	var rows = items.length/cellsPerRow;
	if(tbody.tagName == 'DIV') {		
		var result = "";
		for(var i=0; i < rows; ++i) {
			for(var j=0; j<cellsPerRow; ) {								
				var item = getItem(tmp++);
				if(item.link != '#') {
					if(item.url) {
						result += "<a href='" + item.link + "' title='" + item.title;
						result += "' target='_blank'><img src='" + item.url;
						result += "' style='border: 0' alt='" + item.title + "' /></a>";
					}
					else {
						result += "<b><a href='" + item.link;
						result += "' target='_blank'>" + item.title;
						result += "</a></b>";
					}
					++j;
				}							
				
				if(tmp > items.length)
					break;
			}	
		}
		return result;
	}	
	
	tbody = document.createElement('tbody');
	for(var i=0; i< rows; ++i) {
		var tr = tbody.insertRow(tbody.rows.length);
		for(var j=0; j<cellsPerRow; ) {			
			var item = getItem(tmp++);
			if(item.link != '#') {
				var td = tr.insertCell(tr.cells.length);
				td.setAttribute("width", 100 / cellsPerRow + "%");
				td.setAttribute("height", "24");
				td.bgColor = '#F1F1F1';		
				var result = "";
				if(item.url) {
					result += "<a href='" + item.link + "' title='" + item.title;
					result += "' target='_blank'><img src='" + item.url;
					result += "' style='border: 0' alt='" + item.title + "' /></a>";
				}
				else {
					result += "<b><a href='" + item.link;
					result += "' target='_blank'>" + item.title;
					result += "</a></b>";
				}
				td.innerHTML = result;
				++j;
			}

			if(tmp > items.length)
				return tbody.innerHTML;
		}
	}
	return tbody.innerHTML;
}
// -->
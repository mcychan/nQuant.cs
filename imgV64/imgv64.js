<!--
function init()
{
	cellsPerRow = 5;	
}
function formatDate(pDate)
{
	var months = ['Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
	var dates = pDate.split("/");
	return months[parseInt(dates[1], 10)] + " " + dates[2] + ", " + dates[0];
}

function get_bits_system_architecture()
{
    var _to_check = [] ;
    if (window.navigator.cpuClass) _to_check.push( ( window.navigator.cpuClass + "").toLowerCase()) ;
    if (window.navigator.platform) _to_check.push( ( window.navigator.platform + "").toLowerCase()) ;
    if (navigator.userAgent) _to_check.push( ( navigator.userAgent + "" ).toLowerCase()) ;

    var _64bits_signatures = ["x86_64", "x86-64", "Win64", "x64;", "amd64", "AMD64", "WOW64", "x64_64", "ia64", "sparc64", "ppc64", "IRIX64"] ;
    var _bits = 32, _i, _c;

    for( var _c = 0 ; _c < _to_check.length ; _c++ )
    {
        for( _i = 0 ; _i < _64bits_signatures.length ; _i++ )
        {
            if ( _to_check[_c].indexOf( _64bits_signatures[_i].toLowerCase() ) != -1 )
            {
               _bits = 64;
               return _bits;
            }
        }
    }
    return _bits; 
}

function is_32bits_architecture() { return get_bits_system_architecture() == 32 ? 1 : 0 ; }

function renderVersion(doc)
{
	var current = getNodeNS("", "current", doc, 0);
	var previous = getNodeNS("", "previous", doc, 0);	
	var d1 = formatDate(current.getAttribute("date"));
	var d2 = formatDate(previous.getAttribute("date"));

	var url = getElementTextNS("", "url", doc, 0);
	var url2 = current.firstChild.nodeValue;
	var url64 = getElementTextNS("", "url64", doc, 0);
	var v1 = current.getAttribute("text");
	var v2 = previous.getAttribute("version");
	var curr = current.getAttribute("version");
	
	var flds = ['d1', 'd2', 'v1', 'v2'];
	for(var i=0; i<flds.length; ++i) {
		if(document.getElementById(flds[i]))
			document.getElementById(flds[i]).innerHTML = eval(flds[i]);
	}
	flds = ['url', 'url2'];
	for(var i=0; i<flds.length; ++i) {
		if(document.getElementById(flds[i]))
			document.getElementById(flds[i]).setAttribute("href", eval(flds[i]));
	}
	if(document.getElementById("wincln64"))
		document.getElementById("wincln64").setAttribute("href", url64);

	var pVer = location.search.substring(1, location.search.length);

	if(pVer != '' && curr > pVer)
	{	
		if(confirm('Download new version?'))
			location.href = is_32bits_architecture() ? url : url64;
	}
}

function toggle(tabName) {	 
	var tabs = ['product', 'download', 'support', 'faq'];	
	for(var i=0; i<tabs.length; ++i)		
		document.getElementById(tabs[i]).style.display = (tabName == tabs[i]) ? 'block': 'none';
	
	var tbody = document.getElementById("supporters");
	if(tbody && tbody.innerHTML.trim().length <=0)
		ajax('promotors.xml', renderPage);
}
 //-->
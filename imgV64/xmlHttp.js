<!--
// global flag
var isIE = false;

// retrieve XML document (reusable generic function);
// parameter is URL string (relative or complete) to
// an .xml file whose Content-Type is a valid XML
// type, such as text/xml; XML source must be from
// same domain as HTML file
function loadXMLDoc(url, callback) {
	var xhr = null;
    // branch for native XMLHttpRequest object
    if(window.XMLHttpRequest) {
    	try {
			xhr = new XMLHttpRequest();
        } catch(e) {
			xhr = null;
        }
    // branch for IE/Windows ActiveX version
    } else if(window.ActiveXObject) {
	isIE = true;
       	try {
        	xhr = new ActiveXObject("Msxml2.XMLHTTP");
      	} catch(e) {
        	try {
          		xhr = new ActiveXObject("Microsoft.XMLHTTP");
        	} catch(e) {
          		xhr = null;
        	}
		}
    }
	if(xhr) {
		xhr.onreadystatechange = function() {
			// only if xhr shows "loaded"
			if (xhr.readyState == 4) {
			// only if "OK"
				if (xhr.status == 200 && typeof xhr.callback == "function")
					xhr.callback(xhr.responseXML);
				else
					alert("There was a problem retrieving the XML data:\n" + xhr.statusText);
			}
		};
		xhr.callback = callback;
		xhr.open("GET", url, true);
		xhr.send();
	}
}

// invoked by "Category" select element change;
// loads chosen XML document, clears Topics select
// element, loads new items into Topics select element
function ajax(url, callback) {      
	try { 
		loadXMLDoc(url, callback);	
	}
	catch(e) {
        var msg = (typeof e == "string") ? e : ((e.message) ? e.message : "Unknown Error");
        alert("Unable to get XML data:\n" + msg);
	}
}

// retrieve node of an XML document element, including
// elements using namespaces
function getNodeNS(prefix, local, parentElem, index) {
    if (prefix && isIE) {
        // IE/Windows way of handling namespaces
        result = parentElem.getElementsByTagName(prefix + ":" + local)[index];
    } else {
        // the namespace versions of this method 
        // (getElementsByTagNameNS()) operate
        // differently in Safari and Mozilla, but both
        // return value with just local name, provided 
        // there aren't conflicts with non-namespace element
        // names
        result = parentElem.getElementsByTagName(local)[index];
    }
    return result;
}

// retrieve text of an XML document element, including
// elements using namespaces
function getElementTextNS(prefix, local, parentElem, index) {
    var result = getNodeNS(prefix, local, parentElem, index);
    if (result) {
        // get text, accounting for possible
        // whitespace (carriage return) text nodes 
        if (result.childNodes.length > 1) {
            return result.childNodes[1].nodeValue;
        } else {
            return result.firstChild.nodeValue;    		
        }
    } else {
        return "n/a";
    }
}
//-->
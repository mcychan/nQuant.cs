<!--
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
	var d1 = current.getAttribute("date");
	var d2 = previous.getAttribute("date");

	var url = getElementTextNS("", "url", doc, 0);
	var url2 = current.firstChild.nodeValue;
	var url64 = getElementTextNS("", "url64", doc, 0);
	var v1 = current.getAttribute("text");
	var v2 = previous.getAttribute("version");
	var curr = current.getAttribute("version");
	
	const result = {'d1': d1, 'v1': v1, 'd2': d2, 'v2': v2,
		'url': url, 'url2': url2, 'url64': url64};

	var pVer = location.search.substring(1, location.search.length);

	if(pVer != '' && curr > pVer)
	{	
		if(confirm('Download new version?'))
			location.href = is_32bits_architecture() ? url : url64;
	}
	return result;
}

class Header extends preact.Component {  
	render() {		
		return preact.createElement("div", {className: "container", style: {textAlign: "center"} },
			[
				preact.createElement("h2", {style: {paddingBottom: "1em"}}, preact.createElement("a", {key: "hDownload", name: "download", href: "#", onClick: this.props.onClick}, "ImgV64")),
				preact.createElement("a", {href: "copyright.htm", style: {paddingBottom: "1em"}}, "Copyright \u00a9 2016-2021,"),
				preact.createElement("div", {style: {paddingBottom: "1em"}}, "Miller Cy Chan,"),
				preact.createElement("div", {style: {paddingBottom: "1em"}}, "All rights reserved")
			]
		);
	}
}

class Share extends preact.Component {  
	render() {
		const childrenData = [
			{id: "twitter", attrs: {href: "https://twitter.com/intent/tweet?source=http%3A%2F%2Fimgv64.rf.gd&amp;text=Imgv64%20http%3A%2F%2Fimgv64.rf.gd&amp;via=imgv64"}},
			{id: "google-plus", attrs: {href: "https://plus.google.com/share?url=http%3A%2F%2Fimgv64.rf.gd"}},
			{id: "facebook", attrs: {href: "https://www.facebook.com/sharer/sharer.php?u=http%3A%2F%2Fimgv64.rf.gd&amp;t=Imgv64"}},
			{id: "stumbleupon", attrs: {href: "http://www.stumbleupon.com/submit?url=http%3A%2F%2Fimgv64.rf.gd&amp;title=Imgv64"}},
			{id: "reddit", attrs: {href: "http://www.reddit.com/submit?url=http%3A%2F%2Fimgv64.rf.gd&amp;title=Imgv64"}},
			{id: "linkedin", attrs: {href: "http://www.linkedin.com/shareArticle?mini=true&amp;url=http%3A%2F%2Fimgv64.rf.gd&amp;title=Imgv64&amp;summary=ImgV64%20is%20a%20Windows%20Vista%20or%20above%20graphics%20viewer%20for%20GIF,%20JPG,%20PNG%20and%20other%20formats.&amp;source=http%3A%2F%2Fimgv64.rf.gd"}},
			{id: "envelope", btn: "email", attrs: {href: "mailto:?subject=Imgv64&amp;body=ImgV64%20is%20a%20Windows%20Vista%20or%20above%20graphics%20viewer%20for%20GIF,%20JPG,%20PNG%20and%20other%20formats.%0A%0ARead%20more%20here:%20http%3A%2F%2Fimgv64.rf.gd"}}
		];
		
		return preact.createElement("div", {className: "panel", style: {textAlign: "center", paddingBottom: "1em"}},
			[
				preact.createElement("span", {}, " Share to: "),
				childrenData.map((item, index) => {
					let {btn, id} = item;
					if(!btn)
						btn = id;
					return preact.createElement("a", {key: id, target: "_blank", className: `share-btn ${btn}`, ...item["attrs"]}, 
						preact.createElement("i", {key: `fa-${id}`, className: `fa fa-${id}`, "aria-hidden": true}))
				})
			]
		);
	}
}

class Tabs extends preact.Component {  
	render() {
		const childrenData = [
			{text: "Product Information", attrs: {id: "productInfo", href: "#"}},
			{text: "Download", attrs: {id: "download", href: "#"}},
			{text: "Support", attrs: {id: "support", href: "#"}}
		];
		
		return preact.createElement("div", {className: "panel", style: {paddingBottom: "1em"}},
			[
				childrenData.map((item, index) => {
					return preact.createElement("a", {key: `link${index}`, name: item["id"], onClick: this.props.onClick, ...item["attrs"]},
						preact.createElement("b", {style: {paddingLeft: "2em", paddingRight: "2em"} }, item["text"]))
				})
			]
		);
	}
}

class Promotors extends preact.Component {  
	constructor(props) {
		super(props);
		this.state = { promotors: ''};	
		const tbody = document.getElementById("supporters");
		if(this.state.promotors.trim().length <=0) {
			fetch('promotors.xml') 
			.then(response => response.text()) 
			.then(responseText => new DOMParser().parseFromString(responseText, "text/xml"))
			.then(xmlDoc => {
				const result = renderPage(xmlDoc);
				this.setState({promotors: result}); 
			});
		}
	}
	
	render() {		
		return preact.createElement("div", {id: "supporters", className: "panel", style: {textAlign: "center"},  
			dangerouslySetInnerHTML: { __html:  this.state.promotors} });
	}
}

class Download extends preact.Component {  
	constructor(props) {
		super(props);
		this.state = { d1: '2021/03/20', v1: 1.77, d2: '2020/02/08', v2: 1.76, 
			url: 'ImgV64_17.msi', url2: 'ImgV64_17.msi', url64: 'ImgV64_17_x64.msi'};	
		fetch('history.xml') 
		.then(response => response.text()) 
		.then(responseText => new DOMParser().parseFromString(responseText, "text/xml"))
		.then(xmlDoc => {
			const result = renderVersion(xmlDoc);
			this.setState(result); 
		});					
	}	
	
	render() {	
		const display = this.props.display ? "block" : "none";
		return preact.createElement("div", {id: "download", style: {display: display, paddingBottom: "1em"}},
			[
				preact.createElement("div", {className: "row-padding"}, 
				[
					preact.createElement("div", {className: "half"}, 
					[
						preact.createElement("div", {style: {float: "left", width: "40%", textAlign: "right"}}, 
							preact.createElement("span", {key: "d1"}, formatDate(this.state.d1) )),
						preact.createElement("div", {style: {float: "right", paddingLeft: "2em", width: "50%"}}, 
							preact.createElement("a", {href: this.state.url}, `ImgV64 ${this.state.v1}`))
					]),
					preact.createElement("div", {className: "half"}, 
					[
						preact.createElement("div", {style: {float: "left", width: "40%", textAlign: "right"}}, 
							preact.createElement("span", {key: "d2"}, formatDate(this.state.d2) )),
						preact.createElement("div", {style: {float: "right", paddingLeft: "2em", width: "50%"}}, 
							preact.createElement("a", {href: this.state.url2}, `ImgV64 ${this.state.v2}`))
					])
				]),
				preact.createElement("div", {className: "row-padding"}, 
					preact.createElement("div", {className: "half"}, 
						preact.createElement("div", {style: {clear: "both", textAlign: "center"}}, 
							preact.createElement("a", {href: this.state.url64, id: "ImgV64_17_64"}, "x86 64 bit version Download"))
					)
				)
			],
			preact.createElement(Promotors, {key: "supporters"})
		);
	}
}

class ProductInfo extends preact.Component {  	
	render() {
		const childrenData = [
			[
				'<h3>Introduction</h3><p>Graphics Viewer for Windows Vista or above</p>',
				'<p />ImgV64 is a Windows Vista/Win7/Win8/Win10 graphics viewer for GIF, JPG, PNG, WEBP, HEIC and other formats. It is designed for those who want a solid, easy viewer but who still want powerful features, without the complexity of a complete paint or thumbnailing program.',
				'<p>ImgV64 has many powerful features available, but they won\'t get in your way if you just want to use it as a simple viewer.</p>' +
				'<p>Among the features: You can load an image from a folder and then "slide show" through the rest of the images in the folder using the keyboard left and right arrow keys.</p>' +
				'<p>This works in window or full screen modes.</p>' +
				'<p>You can drag your image from file explorer or browser then drop to ImgV64, or copy url to clipboard then paste onto ImgV64.</p>',
				'<p>Free of any malware, spyware, and viruses. It is a 100% clean and safe tool for you to use.<br/>' +
				'You can add text to an image with the text tool. This is quite versatile, with font, size and color selection, transparent, and 90/180/270 degree rotation of image available.</p>' +
				'<p>As you can see above, many tools, including Undo/Redo, are available on the toolbar. If you want more screen space, you can hide the toolbar and just use menus.</p>'
			],
			[
				'<p>You can choose to open an image at 100% of its original size, or at a size that is adjusted to fill available space on the screen (while maintaining proportions).</p>',
				'<p>You can choose to have up to 4 recently-loaded images in the file menu, for quick re-display. <br />' +
				'If you really want to simply view or edit your photo albums, this is your choice.</p>',
				'<p>Brightness, Contrast, and Gamma adjustments are selectable by slider bar. Other Tune Dialogs include Sharpen, Blur, and color adjustments on the Red/Green/Blue scale.</p>' +
				'<img src="images/ImgV64.jpg" alt="Imgv64 screenshot" class="responsive" />'
			],
			[
				'<p>The Undo of editing also has a Redo option. Both are on the toolbar for easy scrolling through changes.</p>',
				'<p>There is a Set as Wallpaper function in the Image menu to make the current image your desktop wallpaper.</p>',
				'<h3><u>System requirement</u></h3>' +
				'<p>It works on Windows 10, 8 & 8.1, 7, Vista, Server 2008/2012/2016. 64-bit Windows versions are also supported.</p>'
			]
		];
		
		const display = this.props.display ? "block" : "none";
		return preact.createElement("div", {id: "productInfo", style: {display: display, paddingBottom: "1em"}},
			childrenData.map((row, i) => {
				return preact.createElement("div", {className: "row-padding"},
					row.map((text, j) => {
						const className = text.indexOf("<h") >= 0 ? "col" : "third";
						return preact.createElement("div", {key: `info_${j}`, className: className,  
							dangerouslySetInnerHTML: { __html:  text} })
					})
				)
			})
		);
	}
}

class Support extends preact.Component {  	
	render() {
		const childrenData = [
			[
				'<h3>Technical support is free: You can email to <a href="mailto:miller.chan@gmail.com">miller.chan@gmail.com</a></h3>',
				'<p>Please report the problem with description like how to reproduce, screen capture, frequent of occurrance, PC system info.</p>',
				'<p>Currently supported Language: English, traditional' +
				'Chinese, Simplified Chinese, Japanese<br />Thanks for support:',
				'<a title="Japanese.lng" href="mailto:ckh3111@gmail.com">Ricky Chow</a><br />Would you please help to ' +
				'<a href="http://imgv64.rf.gd/english.lng">translate ImgV64</a> by using notepad (Save target as .lng)!<br />',
				'Then <a href="mailto:wincln@usa.com">Email</a> your translation! ' +
				'In return, Imgv64 will link back to your site by the Help -> Online Help menu if user was in your language.</p>'
			]
		];
		
		const display = this.props.display ? "block" : "none";
		return preact.createElement("div", {id: "support", style: {display: display, paddingBottom: "1em"}},
			childrenData.map((row, i) => {
				return preact.createElement("div", {className: "row-padding"},
					row.map((text, j) => {
						return preact.createElement("div", {key: `info_${j}`, className: "third",  
							dangerouslySetInnerHTML: { __html:  text} })
					})
				)
			})
		);
	}
}

class App extends preact.Component {
	constructor(props) {
		super(props);
		this.state = { showDownload: false, showProductInfo: true, showSupport: false};				
	}
	
	componentDidCatch(error, info) {
		console.error(`Error: ${error.message}`);
	}
	
	onClick = (ev) => {
		let name = ev.currentTarget.name;
		if(name == '')
			name = ev.currentTarget.id;
		const tagId = "show" + name[0].toUpperCase() + name.substring(1);
		let displayState = {};
		const showThis = this.state[tagId] ? "none" : "block";
		Object.keys(this.state).map(key => 
			displayState[key] = (key == tagId) ? showThis : this.state[tagId]				
		);
	    this.setState(displayState);
	}
	
	render() {
		return [
			preact.createElement(Header, {key: "header", onClick: (e) => {this.onClick(e)} }),
			preact.createElement(Share, {key: "share"}),
			preact.createElement(Tabs, {key: "tabs", onClick: (e) => {this.onClick(e)} }),
			preact.createElement(Download, {key: "download", display: this.state.showDownload}),
			preact.createElement(ProductInfo, {key: "productInfo", display: this.state.showProductInfo}),
			preact.createElement(Support, {key: "support", display: this.state.showSupport})
		];
	}
}

// render
preact.render(preact.createElement(App, {}), document.querySelector('#imgv64'));
 //-->
var NAVTREE =
[
  [ "Matrix", "index.html", [
    [ "Data Structures", "annotated.html", [
      [ "cs_dmperm_results", "structcs__dmperm__results.html", null ],
      [ "cs_numeric", "structcs__numeric.html", null ],
      [ "cs_sparse", "structcs__sparse.html", null ],
      [ "cs_symbolic", "structcs__symbolic.html", null ]
    ] ],
    [ "Data Structure Index", "classes.html", null ],
    [ "Data Fields", "functions.html", null ],
    [ "File List", "files.html", [
      [ "abIndex.c", "abIndex_8c.html", null ],
      [ "abIndex.h", "abIndex_8h.html", null ],
      [ "chm_common.c", "chm__common_8c.html", null ],
      [ "chm_common.h", "chm__common_8h.html", null ],
      [ "CHMfactor.c", "CHMfactor_8c.html", null ],
      [ "CHMfactor.h", "CHMfactor_8h.html", null ],
      [ "cs.c", "cs_8c.html", null ],
      [ "cs.h", "cs_8h.html", null ],
      [ "cs_utils.c", "cs__utils_8c.html", null ],
      [ "cs_utils.h", "cs__utils_8h.html", null ],
      [ "Csparse.c", "Csparse_8c.html", null ],
      [ "Csparse.h", "Csparse_8h.html", null ],
      [ "dense.c", "dense_8c.html", null ],
      [ "dense.h", "dense_8h.html", null ],
      [ "dgCMatrix.c", "dgCMatrix_8c.html", null ],
      [ "dgCMatrix.h", "dgCMatrix_8h.html", null ],
      [ "dgeMatrix.c", "dgeMatrix_8c.html", null ],
      [ "dgeMatrix.h", "dgeMatrix_8h.html", null ],
      [ "dgTMatrix.c", "dgTMatrix_8c.html", null ],
      [ "dgTMatrix.h", "dgTMatrix_8h.html", null ],
      [ "dpoMatrix.c", "dpoMatrix_8c.html", null ],
      [ "dpoMatrix.h", "dpoMatrix_8h.html", null ],
      [ "dppMatrix.c", "dppMatrix_8c.html", null ],
      [ "dppMatrix.h", "dppMatrix_8h.html", null ],
      [ "dsCMatrix.c", "dsCMatrix_8c.html", null ],
      [ "dsCMatrix.h", "dsCMatrix_8h.html", null ],
      [ "dspMatrix.c", "dspMatrix_8c.html", null ],
      [ "dspMatrix.h", "dspMatrix_8h.html", null ],
      [ "dsyMatrix.c", "dsyMatrix_8c.html", null ],
      [ "dsyMatrix.h", "dsyMatrix_8h.html", null ],
      [ "dtCMatrix.c", "dtCMatrix_8c.html", null ],
      [ "dtCMatrix.h", "dtCMatrix_8h.html", null ],
      [ "dtpMatrix.c", "dtpMatrix_8c.html", null ],
      [ "dtpMatrix.h", "dtpMatrix_8h.html", null ],
      [ "dtrMatrix.c", "dtrMatrix_8c.html", null ],
      [ "dtrMatrix.h", "dtrMatrix_8h.html", null ],
      [ "dtTMatrix.c", "dtTMatrix_8c.html", null ],
      [ "dtTMatrix.h", "dtTMatrix_8h.html", null ],
      [ "factorizations.c", "factorizations_8c.html", null ],
      [ "factorizations.h", "factorizations_8h.html", null ],
      [ "init.c", "init_8c.html", null ],
      [ "ldense.c", "ldense_8c.html", null ],
      [ "ldense.h", "ldense_8h.html", null ],
      [ "lgCMatrix.c", "lgCMatrix_8c.html", null ],
      [ "lgCMatrix.h", "lgCMatrix_8h.html", null ],
      [ "Mutils.c", "Mutils_8c.html", null ],
      [ "Mutils.h", "Mutils_8h.html", null ],
      [ "sparseQR.c", "sparseQR_8c.html", null ],
      [ "sparseQR.h", "sparseQR_8h.html", null ],
      [ "Syms.h", "Syms_8h.html", null ],
      [ "t_Csparse_subassign.c", "t__Csparse__subassign_8c.html", null ],
      [ "t_gCMatrix_colSums.c", "t__gCMatrix__colSums_8c.html", null ],
      [ "t_Matrix_rle.c", "t__Matrix__rle_8c.html", null ],
      [ "t_sparseVector.c", "t__sparseVector_8c.html", null ],
      [ "TMatrix_as.c", "TMatrix__as_8c.html", null ],
      [ "TMatrix_as.h", "TMatrix__as_8h.html", null ],
      [ "Tsparse.c", "Tsparse_8c.html", null ],
      [ "Tsparse.h", "Tsparse_8h.html", null ]
    ] ],
    [ "Globals", "globals.html", null ]
  ] ]
];

function createIndent(o,domNode,node,level)
{
  if (node.parentNode && node.parentNode.parentNode)
  {
    createIndent(o,domNode,node.parentNode,level+1);
  }
  var imgNode = document.createElement("img");
  if (level==0 && node.childrenData)
  {
    node.plus_img = imgNode;
    node.expandToggle = document.createElement("a");
    node.expandToggle.href = "javascript:void(0)";
    node.expandToggle.onclick = function() 
    {
      if (node.expanded) 
      {
        $(node.getChildrenUL()).slideUp("fast");
        if (node.isLast)
        {
          node.plus_img.src = node.relpath+"ftv2plastnode.png";
        }
        else
        {
          node.plus_img.src = node.relpath+"ftv2pnode.png";
        }
        node.expanded = false;
      } 
      else 
      {
        expandNode(o, node, false);
      }
    }
    node.expandToggle.appendChild(imgNode);
    domNode.appendChild(node.expandToggle);
  }
  else
  {
    domNode.appendChild(imgNode);
  }
  if (level==0)
  {
    if (node.isLast)
    {
      if (node.childrenData)
      {
        imgNode.src = node.relpath+"ftv2plastnode.png";
      }
      else
      {
        imgNode.src = node.relpath+"ftv2lastnode.png";
        domNode.appendChild(imgNode);
      }
    }
    else
    {
      if (node.childrenData)
      {
        imgNode.src = node.relpath+"ftv2pnode.png";
      }
      else
      {
        imgNode.src = node.relpath+"ftv2node.png";
        domNode.appendChild(imgNode);
      }
    }
  }
  else
  {
    if (node.isLast)
    {
      imgNode.src = node.relpath+"ftv2blank.png";
    }
    else
    {
      imgNode.src = node.relpath+"ftv2vertline.png";
    }
  }
  imgNode.border = "0";
}

function newNode(o, po, text, link, childrenData, lastNode)
{
  var node = new Object();
  node.children = Array();
  node.childrenData = childrenData;
  node.depth = po.depth + 1;
  node.relpath = po.relpath;
  node.isLast = lastNode;

  node.li = document.createElement("li");
  po.getChildrenUL().appendChild(node.li);
  node.parentNode = po;

  node.itemDiv = document.createElement("div");
  node.itemDiv.className = "item";

  node.labelSpan = document.createElement("span");
  node.labelSpan.className = "label";

  createIndent(o,node.itemDiv,node,0);
  node.itemDiv.appendChild(node.labelSpan);
  node.li.appendChild(node.itemDiv);

  var a = document.createElement("a");
  node.labelSpan.appendChild(a);
  node.label = document.createTextNode(text);
  a.appendChild(node.label);
  if (link) 
  {
    a.href = node.relpath+link;
  } 
  else 
  {
    if (childrenData != null) 
    {
      a.className = "nolink";
      a.href = "javascript:void(0)";
      a.onclick = node.expandToggle.onclick;
      node.expanded = false;
    }
  }

  node.childrenUL = null;
  node.getChildrenUL = function() 
  {
    if (!node.childrenUL) 
    {
      node.childrenUL = document.createElement("ul");
      node.childrenUL.className = "children_ul";
      node.childrenUL.style.display = "none";
      node.li.appendChild(node.childrenUL);
    }
    return node.childrenUL;
  };

  return node;
}

function showRoot()
{
  var headerHeight = $("#top").height();
  var footerHeight = $("#nav-path").height();
  var windowHeight = $(window).height() - headerHeight - footerHeight;
  navtree.scrollTo('#selected',0,{offset:-windowHeight/2});
}

function expandNode(o, node, imm)
{
  if (node.childrenData && !node.expanded) 
  {
    if (!node.childrenVisited) 
    {
      getNode(o, node);
    }
    if (imm)
    {
      $(node.getChildrenUL()).show();
    } 
    else 
    {
      $(node.getChildrenUL()).slideDown("fast",showRoot);
    }
    if (node.isLast)
    {
      node.plus_img.src = node.relpath+"ftv2mlastnode.png";
    }
    else
    {
      node.plus_img.src = node.relpath+"ftv2mnode.png";
    }
    node.expanded = true;
  }
}

function getNode(o, po)
{
  po.childrenVisited = true;
  var l = po.childrenData.length-1;
  for (var i in po.childrenData) 
  {
    var nodeData = po.childrenData[i];
    po.children[i] = newNode(o, po, nodeData[0], nodeData[1], nodeData[2],
        i==l);
  }
}

function findNavTreePage(url, data)
{
  var nodes = data;
  var result = null;
  for (var i in nodes) 
  {
    var d = nodes[i];
    if (d[1] == url) 
    {
      return new Array(i);
    }
    else if (d[2] != null) // array of children
    {
      result = findNavTreePage(url, d[2]);
      if (result != null) 
      {
        return (new Array(i).concat(result));
      }
    }
  }
  return null;
}

function initNavTree(toroot,relpath)
{
  var o = new Object();
  o.toroot = toroot;
  o.node = new Object();
  o.node.li = document.getElementById("nav-tree-contents");
  o.node.childrenData = NAVTREE;
  o.node.children = new Array();
  o.node.childrenUL = document.createElement("ul");
  o.node.getChildrenUL = function() { return o.node.childrenUL; };
  o.node.li.appendChild(o.node.childrenUL);
  o.node.depth = 0;
  o.node.relpath = relpath;

  getNode(o, o.node);

  o.breadcrumbs = findNavTreePage(toroot, NAVTREE);
  if (o.breadcrumbs == null)
  {
    o.breadcrumbs = findNavTreePage("index.html",NAVTREE);
  }
  if (o.breadcrumbs != null && o.breadcrumbs.length>0)
  {
    var p = o.node;
    for (var i in o.breadcrumbs) 
    {
      var j = o.breadcrumbs[i];
      p = p.children[j];
      expandNode(o,p,true);
    }
    p.itemDiv.className = p.itemDiv.className + " selected";
    p.itemDiv.id = "selected";
    $(window).load(showRoot);
  }
}


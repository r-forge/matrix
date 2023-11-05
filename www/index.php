
<!-- This is == http://matrix.r-forge.r-project.org/ -->
<!--            ==================================== -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><The Matrix Project></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h1>The Matrix Project</h1>

<p>
This is the front page of <a href="https://r-forge.r-project.org/projects/matrix/">the Matrix project</a> on R-Forge, which is the home of R packages <a href="https://cran.r-project.org/package=Matrix"><b>Matrix</b></a> and <a href="https://cran.r-project.org/package=MatrixModels"><b>MatrixModels</b></a>.  These provide S4 classes and methods for representing, operating on, and modelling with sparse and dense matrices.
</p>

<h2>Installation</h2>

<p>
The latest release and development versions of <b>Matrix</b> can be installed from CRAN and R-Forge, respectively, with:
</p>

<pre>
> install.packages("Matrix")
> install.packages("Matrix", repos = "http://R-Forge.R-project.org")
</pre>

<p>
Older release versions are preserved in the CRAN <a href="https://cran.r-project.org/src/contrib/Archive/Matrix/">archive</a>.
</p>

<p>
Since R version 2.9.0, <b>Matrix</b> has belonged to a list of "recommended" packages bundled in standard installations of R, so that <code>library(Matrix)</code> works out of the box.  Still, one often installs <b>Matrix</b> from CRAN into a user library if the version in the default library is out of date.
</p>

<p>
Note that <em>source</em> installation (as opposed to <em>binary</em> installation) requires compilers and additional tools; see the <a href="https://stat.ethz.ch/R-manual/R-patched/doc/manual/R-admin.html#Installing-packages"><i>R Installation and Administration</i></a> manual for platform-specific details.
</p>

<h2>Documentation</h2>

<p>
Indices of available help topics, data sets, and vignettes can be accessed with:
</p>

<pre>
> help(package = "Matrix")
> data(package = "Matrix")
> vignette(package = "Matrix")
</pre>

<p>
Passing <code>help_type = "html"</code> to <code>help</code> renders a hyperlinked version of the index in your browser.  Because <b>Matrix</b> is "recommended", its index is also hosted <a href="https://stat.ethz.ch/R-manual/R-patched/library/Matrix/html/00Index.html">here</a>.  (The version is that which is bundled in the latest release of R, which will be recent but need not be current.)
</p>

<p>
Doxygen-generated documentation for the C sources of <b>Matrix</b> is hosted <a href="./doxygen">here</a>.  Do nudge <code>maintainer("Matrix")</code> if it seems to be severely out of date.
</p>

Slides for past talks about the Matrix project are hosted <a href="./slides">here</a>.

<h2>Accessing the source code</h2>

<p>
The Matrix project is maintained in a Subversion ("SVN") repository.  If you have installed Subversion, then you can check out the latest development version of <b>Matrix</b> with:
</p>

<pre>
$ svn checkout [-r REV] svn://svn.r-forge.r-project.org/svnroot/matrix/pkg/Matrix [DEST]
</pre>

<p>
where <code>REV</code> and <code>DEST</code> optionally specify a revision number and destination path.  Without Subversion, one can browse the repository using R-Forge's <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=matrix">ViewVC interface</a> or download a tarball produced by <code>R CMD build</code> from CRAN or R-Forge:
</p>

<pre>
> download.packages("Matrix", ".", type = "source")
> download.packages("Matrix", ".", type = "source", repos = "http://R-Forge.R-project.org")
</pre>

<h2>Bug reports and feature requests</h2>

The Matrix project uses <a href="https://r-forge.r-project.org/tracker/?group_id=61">trackers</a> provided by R-Forge to manage bug reports and feature requests.  An R-Forge user account is required to post.  Once logged in, select an appropriate tracker and click "Submit New".  Before posting, <em>do</em> verify that your issue exists under the latest development version of the relevant package.

</body>
</html>

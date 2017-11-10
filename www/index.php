
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
	<title><?php echo $group_name; ?></title>
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

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<p>
<t>Matrix</t> has been a "recommended" R package since R version 2.9.0 (October 2008). This means that it is now part of every regular installation of R.
</p>

<p>
The <t>Matrix</t> package comes with several vignettes describing its functionality.
In addition, we provide <strong><a href="./slides">slides</a></strong> of presentations or talks about the subject.
</p>

<p> <strong>Doxygen</strong> documentation of the underlying C functions is <a href="./doxygen/"><strong>here</strong></a>. </p>

<p> <strong>Do read</string> on <strong><a href="./programming/">compatibility</a></strong> changes  in connection with <t>Matrix</t>.
</p>


<p> The <strong>project summary page</strong> is <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>.

    From there you can browse the source via viewVC, for the Matrix package more directly
    <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/Matrix/?root=matrix"><strong>here</strong></a>,

    or look at the build and `R CMD check` status and also download the result,
    <a href="https://r-forge.r-project.org/R/?group_id=61" <strong>here</strong></a>.
</p>

</body>
</html>

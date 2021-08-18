cola Report for Consensus Partitioning
==================

**Date**: 2021-01-11 00:37:22 CST, **cola version**: 1.9.0.1013

----------------------------------------------------------------

<style type='text/css'>

body, td, th {
   font-family: Arial,Helvetica,sans-serif;
   background-color: white;
   font-size: 13px;
  max-width: 800px;
  margin: auto;
  margin-left:210px;
  padding: 0px 10px 0px 10px;
  border-left: 1px solid #EEEEEE;
  line-height: 150%;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, 

monospace;
}

h1 {
   font-size:2.2em;
}

h2 {
   font-size:1.8em;
}

h3 {
   font-size:1.4em;
}

h4 {
   font-size:1.0em;
}

h5 {
   font-size:0.9em;
}

h6 {
   font-size:0.8em;
}

a {
  text-decoration: none;
  color: #0366d6;
}

a:hover {
  text-decoration: underline;
}

a:visited {
   color: #0366d6;
}

pre, img {
  max-width: 100%;
}
pre {
  overflow-x: auto;
}
pre code {
   display: block; padding: 0.5em;
}

code {
  font-size: 92%;
  border: 1px solid #ccc;
}

code[class] {
  background-color: #F8F8F8;
}

table, td, th {
  border: 1px solid #ccc;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>




## Summary



First the variable is renamed to `res_list`.


```r
res_list = final
```



All available functions which can be applied to this `res_list` object:


```r
res_list
```

```
#> A 'ConsensusPartitionList' object with 1 methods.
#>   On a matrix with 13 rows and 2737 columns.
#>   Top rows are extracted by 'ATC' methods.
#>   Subgroups are detected by 'skmeans' method.
#>   Number of partitions are tried for k = 2, 3, 4, 5, 6, 7, 8, 9, 10.
#>   Performed in total 450 partitions by row resampling.
#> 
#> Following methods can be applied to this 'ConsensusPartitionList' object:
#>  [1] "cola_report"           "collect_classes"       "collect_plots"         "collect_stats"        
#>  [5] "colnames"              "functional_enrichment" "get_anno_col"          "get_anno"             
#>  [9] "get_classes"           "get_matrix"            "get_membership"        "get_stats"            
#> [13] "is_best_k"             "is_stable_k"           "ncol"                  "nrow"                 
#> [17] "rownames"              "show"                  "suggest_best_k"        "test_to_known_factors"
#> [21] "top_rows_heatmap"      "top_rows_overlap"     
#> 
#> You can get result for a single method by, e.g. object["ATC", "skmeans"] or object["ATC:skmeans"]
```

The call of `run_all_consensus_partition_methods()` was:


```
#> run_all_consensus_partition_methods(data = mat_adj, top_value_method = "ATC", partition_method = "skmeans", 
#>     max_k = 10, top_n = 13, mc.cores = 8)
```

Dimension of the input matrix:


```r
mat = get_matrix(res_list)
dim(mat)
```

```
#> [1]   13 2737
```

### Density distribution

The density distribution for each sample is visualized as in one column in the
following heatmap. The clustering is based on the distance which is the
Kolmogorov-Smirnov statistic between two distributions.




```r
library(ComplexHeatmap)
densityHeatmap(mat, ylab = "value", cluster_columns = TRUE, show_column_names = FALSE,
    mc.cores = 1)
```





### Suggest the best k



Folowing table shows the best `k` (number of partitions) for each combination
of top-value methods and partitioning methods. Clicking on the method name in
the table goes to the corresponding section for a single combination of methods.

[The cola vignette](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_13)
explains the definition of the metrics used for determining the best
number of partitions.


```r
suggest_best_k(res_list)
```


|                            | The best k| 1-PAC| Mean silhouette| Concordance|   |
|:---------------------------|----------:|-----:|---------------:|-----------:|:--|
|[ATC:skmeans](#ATC-skmeans) |          2| 0.679|           0.858|       0.937|   |

\*\*: 1-PAC > 0.95, \*: 1-PAC > 0.9




### CDF of consensus matrices

Cumulative distribution function curves of consensus matrix for all methods.




```r
collect_plots(res_list, fun = plot_ecdf)
```



### Consensus heatmap

Consensus heatmaps for all methods. ([What is a consensus heatmap?](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_9))


<style type='text/css'>



.ui-helper-hidden {
	display: none;
}
.ui-helper-hidden-accessible {
	border: 0;
	clip: rect(0 0 0 0);
	height: 1px;
	margin: -1px;
	overflow: hidden;
	padding: 0;
	position: absolute;
	width: 1px;
}
.ui-helper-reset {
	margin: 0;
	padding: 0;
	border: 0;
	outline: 0;
	line-height: 1.3;
	text-decoration: none;
	font-size: 100%;
	list-style: none;
}
.ui-helper-clearfix:before,
.ui-helper-clearfix:after {
	content: "";
	display: table;
	border-collapse: collapse;
}
.ui-helper-clearfix:after {
	clear: both;
}
.ui-helper-zfix {
	width: 100%;
	height: 100%;
	top: 0;
	left: 0;
	position: absolute;
	opacity: 0;
	filter:Alpha(Opacity=0); 
}

.ui-front {
	z-index: 100;
}



.ui-state-disabled {
	cursor: default !important;
	pointer-events: none;
}



.ui-icon {
	display: inline-block;
	vertical-align: middle;
	margin-top: -.25em;
	position: relative;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
}

.ui-widget-icon-block {
	left: 50%;
	margin-left: -8px;
	display: block;
}




.ui-widget-overlay {
	position: fixed;
	top: 0;
	left: 0;
	width: 100%;
	height: 100%;
}
.ui-accordion .ui-accordion-header {
	display: block;
	cursor: pointer;
	position: relative;
	margin: 2px 0 0 0;
	padding: .5em .5em .5em .7em;
	font-size: 100%;
}
.ui-accordion .ui-accordion-content {
	padding: 1em 2.2em;
	border-top: 0;
	overflow: auto;
}
.ui-autocomplete {
	position: absolute;
	top: 0;
	left: 0;
	cursor: default;
}
.ui-menu {
	list-style: none;
	padding: 0;
	margin: 0;
	display: block;
	outline: 0;
}
.ui-menu .ui-menu {
	position: absolute;
}
.ui-menu .ui-menu-item {
	margin: 0;
	cursor: pointer;
	
	list-style-image: url("data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7");
}
.ui-menu .ui-menu-item-wrapper {
	position: relative;
	padding: 3px 1em 3px .4em;
}
.ui-menu .ui-menu-divider {
	margin: 5px 0;
	height: 0;
	font-size: 0;
	line-height: 0;
	border-width: 1px 0 0 0;
}
.ui-menu .ui-state-focus,
.ui-menu .ui-state-active {
	margin: -1px;
}


.ui-menu-icons {
	position: relative;
}
.ui-menu-icons .ui-menu-item-wrapper {
	padding-left: 2em;
}


.ui-menu .ui-icon {
	position: absolute;
	top: 0;
	bottom: 0;
	left: .2em;
	margin: auto 0;
}


.ui-menu .ui-menu-icon {
	left: auto;
	right: 0;
}
.ui-button {
	padding: .4em 1em;
	display: inline-block;
	position: relative;
	line-height: normal;
	margin-right: .1em;
	cursor: pointer;
	vertical-align: middle;
	text-align: center;
	-webkit-user-select: none;
	-moz-user-select: none;
	-ms-user-select: none;
	user-select: none;

	
	overflow: visible;
}

.ui-button,
.ui-button:link,
.ui-button:visited,
.ui-button:hover,
.ui-button:active {
	text-decoration: none;
}


.ui-button-icon-only {
	width: 2em;
	box-sizing: border-box;
	text-indent: -9999px;
	white-space: nowrap;
}


input.ui-button.ui-button-icon-only {
	text-indent: 0;
}


.ui-button-icon-only .ui-icon {
	position: absolute;
	top: 50%;
	left: 50%;
	margin-top: -8px;
	margin-left: -8px;
}

.ui-button.ui-icon-notext .ui-icon {
	padding: 0;
	width: 2.1em;
	height: 2.1em;
	text-indent: -9999px;
	white-space: nowrap;

}

input.ui-button.ui-icon-notext .ui-icon {
	width: auto;
	height: auto;
	text-indent: 0;
	white-space: normal;
	padding: .4em 1em;
}



input.ui-button::-moz-focus-inner,
button.ui-button::-moz-focus-inner {
	border: 0;
	padding: 0;
}
.ui-controlgroup {
	vertical-align: middle;
	display: inline-block;
}
.ui-controlgroup > .ui-controlgroup-item {
	float: left;
	margin-left: 0;
	margin-right: 0;
}
.ui-controlgroup > .ui-controlgroup-item:focus,
.ui-controlgroup > .ui-controlgroup-item.ui-visual-focus {
	z-index: 9999;
}
.ui-controlgroup-vertical > .ui-controlgroup-item {
	display: block;
	float: none;
	width: 100%;
	margin-top: 0;
	margin-bottom: 0;
	text-align: left;
}
.ui-controlgroup-vertical .ui-controlgroup-item {
	box-sizing: border-box;
}
.ui-controlgroup .ui-controlgroup-label {
	padding: .4em 1em;
}
.ui-controlgroup .ui-controlgroup-label span {
	font-size: 80%;
}
.ui-controlgroup-horizontal .ui-controlgroup-label + .ui-controlgroup-item {
	border-left: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label + .ui-controlgroup-item {
	border-top: none;
}
.ui-controlgroup-horizontal .ui-controlgroup-label.ui-widget-content {
	border-right: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label.ui-widget-content {
	border-bottom: none;
}


.ui-controlgroup-vertical .ui-spinner-input {

	
	width: 75%;
	width: calc( 100% - 2.4em );
}
.ui-controlgroup-vertical .ui-spinner .ui-spinner-up {
	border-top-style: solid;
}

.ui-checkboxradio-label .ui-icon-background {
	box-shadow: inset 1px 1px 1px #ccc;
	border-radius: .12em;
	border: none;
}
.ui-checkboxradio-radio-label .ui-icon-background {
	width: 16px;
	height: 16px;
	border-radius: 1em;
	overflow: visible;
	border: none;
}
.ui-checkboxradio-radio-label.ui-checkboxradio-checked .ui-icon,
.ui-checkboxradio-radio-label.ui-checkboxradio-checked:hover .ui-icon {
	background-image: none;
	width: 8px;
	height: 8px;
	border-width: 4px;
	border-style: solid;
}
.ui-checkboxradio-disabled {
	pointer-events: none;
}
.ui-datepicker {
	width: 17em;
	padding: .2em .2em 0;
	display: none;
}
.ui-datepicker .ui-datepicker-header {
	position: relative;
	padding: .2em 0;
}
.ui-datepicker .ui-datepicker-prev,
.ui-datepicker .ui-datepicker-next {
	position: absolute;
	top: 2px;
	width: 1.8em;
	height: 1.8em;
}
.ui-datepicker .ui-datepicker-prev-hover,
.ui-datepicker .ui-datepicker-next-hover {
	top: 1px;
}
.ui-datepicker .ui-datepicker-prev {
	left: 2px;
}
.ui-datepicker .ui-datepicker-next {
	right: 2px;
}
.ui-datepicker .ui-datepicker-prev-hover {
	left: 1px;
}
.ui-datepicker .ui-datepicker-next-hover {
	right: 1px;
}
.ui-datepicker .ui-datepicker-prev span,
.ui-datepicker .ui-datepicker-next span {
	display: block;
	position: absolute;
	left: 50%;
	margin-left: -8px;
	top: 50%;
	margin-top: -8px;
}
.ui-datepicker .ui-datepicker-title {
	margin: 0 2.3em;
	line-height: 1.8em;
	text-align: center;
}
.ui-datepicker .ui-datepicker-title select {
	font-size: 1em;
	margin: 1px 0;
}
.ui-datepicker select.ui-datepicker-month,
.ui-datepicker select.ui-datepicker-year {
	width: 45%;
}
.ui-datepicker table {
	width: 100%;
	font-size: .9em;
	border-collapse: collapse;
	margin: 0 0 .4em;
}
.ui-datepicker th {
	padding: .7em .3em;
	text-align: center;
	font-weight: bold;
	border: 0;
}
.ui-datepicker td {
	border: 0;
	padding: 1px;
}
.ui-datepicker td span,
.ui-datepicker td a {
	display: block;
	padding: .2em;
	text-align: right;
	text-decoration: none;
}
.ui-datepicker .ui-datepicker-buttonpane {
	background-image: none;
	margin: .7em 0 0 0;
	padding: 0 .2em;
	border-left: 0;
	border-right: 0;
	border-bottom: 0;
}
.ui-datepicker .ui-datepicker-buttonpane button {
	float: right;
	margin: .5em .2em .4em;
	cursor: pointer;
	padding: .2em .6em .3em .6em;
	width: auto;
	overflow: visible;
}
.ui-datepicker .ui-datepicker-buttonpane button.ui-datepicker-current {
	float: left;
}


.ui-datepicker.ui-datepicker-multi {
	width: auto;
}
.ui-datepicker-multi .ui-datepicker-group {
	float: left;
}
.ui-datepicker-multi .ui-datepicker-group table {
	width: 95%;
	margin: 0 auto .4em;
}
.ui-datepicker-multi-2 .ui-datepicker-group {
	width: 50%;
}
.ui-datepicker-multi-3 .ui-datepicker-group {
	width: 33.3%;
}
.ui-datepicker-multi-4 .ui-datepicker-group {
	width: 25%;
}
.ui-datepicker-multi .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-multi .ui-datepicker-group-middle .ui-datepicker-header {
	border-left-width: 0;
}
.ui-datepicker-multi .ui-datepicker-buttonpane {
	clear: left;
}
.ui-datepicker-row-break {
	clear: both;
	width: 100%;
	font-size: 0;
}


.ui-datepicker-rtl {
	direction: rtl;
}
.ui-datepicker-rtl .ui-datepicker-prev {
	right: 2px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next {
	left: 2px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-prev:hover {
	right: 1px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next:hover {
	left: 1px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane {
	clear: right;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button {
	float: left;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button.ui-datepicker-current,
.ui-datepicker-rtl .ui-datepicker-group {
	float: right;
}
.ui-datepicker-rtl .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-rtl .ui-datepicker-group-middle .ui-datepicker-header {
	border-right-width: 0;
	border-left-width: 1px;
}


.ui-datepicker .ui-icon {
	display: block;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
	left: .5em;
	top: .3em;
}
.ui-dialog {
	position: absolute;
	top: 0;
	left: 0;
	padding: .2em;
	outline: 0;
}
.ui-dialog .ui-dialog-titlebar {
	padding: .4em 1em;
	position: relative;
}
.ui-dialog .ui-dialog-title {
	float: left;
	margin: .1em 0;
	white-space: nowrap;
	width: 90%;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-dialog .ui-dialog-titlebar-close {
	position: absolute;
	right: .3em;
	top: 50%;
	width: 20px;
	margin: -10px 0 0 0;
	padding: 1px;
	height: 20px;
}
.ui-dialog .ui-dialog-content {
	position: relative;
	border: 0;
	padding: .5em 1em;
	background: none;
	overflow: auto;
}
.ui-dialog .ui-dialog-buttonpane {
	text-align: left;
	border-width: 1px 0 0 0;
	background-image: none;
	margin-top: .5em;
	padding: .3em 1em .5em .4em;
}
.ui-dialog .ui-dialog-buttonpane .ui-dialog-buttonset {
	float: right;
}
.ui-dialog .ui-dialog-buttonpane button {
	margin: .5em .4em .5em 0;
	cursor: pointer;
}
.ui-dialog .ui-resizable-n {
	height: 2px;
	top: 0;
}
.ui-dialog .ui-resizable-e {
	width: 2px;
	right: 0;
}
.ui-dialog .ui-resizable-s {
	height: 2px;
	bottom: 0;
}
.ui-dialog .ui-resizable-w {
	width: 2px;
	left: 0;
}
.ui-dialog .ui-resizable-se,
.ui-dialog .ui-resizable-sw,
.ui-dialog .ui-resizable-ne,
.ui-dialog .ui-resizable-nw {
	width: 7px;
	height: 7px;
}
.ui-dialog .ui-resizable-se {
	right: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-sw {
	left: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-ne {
	right: 0;
	top: 0;
}
.ui-dialog .ui-resizable-nw {
	left: 0;
	top: 0;
}
.ui-draggable .ui-dialog-titlebar {
	cursor: move;
}
.ui-draggable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable {
	position: relative;
}
.ui-resizable-handle {
	position: absolute;
	font-size: 0.1px;
	display: block;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable-disabled .ui-resizable-handle,
.ui-resizable-autohide .ui-resizable-handle {
	display: none;
}
.ui-resizable-n {
	cursor: n-resize;
	height: 7px;
	width: 100%;
	top: -5px;
	left: 0;
}
.ui-resizable-s {
	cursor: s-resize;
	height: 7px;
	width: 100%;
	bottom: -5px;
	left: 0;
}
.ui-resizable-e {
	cursor: e-resize;
	width: 7px;
	right: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-w {
	cursor: w-resize;
	width: 7px;
	left: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-se {
	cursor: se-resize;
	width: 12px;
	height: 12px;
	right: 1px;
	bottom: 1px;
}
.ui-resizable-sw {
	cursor: sw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	bottom: -5px;
}
.ui-resizable-nw {
	cursor: nw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	top: -5px;
}
.ui-resizable-ne {
	cursor: ne-resize;
	width: 9px;
	height: 9px;
	right: -5px;
	top: -5px;
}
.ui-progressbar {
	height: 2em;
	text-align: left;
	overflow: hidden;
}
.ui-progressbar .ui-progressbar-value {
	margin: -1px;
	height: 100%;
}
.ui-progressbar .ui-progressbar-overlay {
	background: url("data:image/gif;base64,R0lGODlhKAAoAIABAAAAAP///yH/C05FVFNDQVBFMi4wAwEAAAAh+QQJAQABACwAAAAAKAAoAAACkYwNqXrdC52DS06a7MFZI+4FHBCKoDeWKXqymPqGqxvJrXZbMx7Ttc+w9XgU2FB3lOyQRWET2IFGiU9m1frDVpxZZc6bfHwv4c1YXP6k1Vdy292Fb6UkuvFtXpvWSzA+HycXJHUXiGYIiMg2R6W459gnWGfHNdjIqDWVqemH2ekpObkpOlppWUqZiqr6edqqWQAAIfkECQEAAQAsAAAAACgAKAAAApSMgZnGfaqcg1E2uuzDmmHUBR8Qil95hiPKqWn3aqtLsS18y7G1SzNeowWBENtQd+T1JktP05nzPTdJZlR6vUxNWWjV+vUWhWNkWFwxl9VpZRedYcflIOLafaa28XdsH/ynlcc1uPVDZxQIR0K25+cICCmoqCe5mGhZOfeYSUh5yJcJyrkZWWpaR8doJ2o4NYq62lAAACH5BAkBAAEALAAAAAAoACgAAAKVDI4Yy22ZnINRNqosw0Bv7i1gyHUkFj7oSaWlu3ovC8GxNso5fluz3qLVhBVeT/Lz7ZTHyxL5dDalQWPVOsQWtRnuwXaFTj9jVVh8pma9JjZ4zYSj5ZOyma7uuolffh+IR5aW97cHuBUXKGKXlKjn+DiHWMcYJah4N0lYCMlJOXipGRr5qdgoSTrqWSq6WFl2ypoaUAAAIfkECQEAAQAsAAAAACgAKAAAApaEb6HLgd/iO7FNWtcFWe+ufODGjRfoiJ2akShbueb0wtI50zm02pbvwfWEMWBQ1zKGlLIhskiEPm9R6vRXxV4ZzWT2yHOGpWMyorblKlNp8HmHEb/lCXjcW7bmtXP8Xt229OVWR1fod2eWqNfHuMjXCPkIGNileOiImVmCOEmoSfn3yXlJWmoHGhqp6ilYuWYpmTqKUgAAIfkECQEAAQAsAAAAACgAKAAAApiEH6kb58biQ3FNWtMFWW3eNVcojuFGfqnZqSebuS06w5V80/X02pKe8zFwP6EFWOT1lDFk8rGERh1TTNOocQ61Hm4Xm2VexUHpzjymViHrFbiELsefVrn6XKfnt2Q9G/+Xdie499XHd2g4h7ioOGhXGJboGAnXSBnoBwKYyfioubZJ2Hn0RuRZaflZOil56Zp6iioKSXpUAAAh+QQJAQABACwAAAAAKAAoAAACkoQRqRvnxuI7kU1a1UU5bd5tnSeOZXhmn5lWK3qNTWvRdQxP8qvaC+/yaYQzXO7BMvaUEmJRd3TsiMAgswmNYrSgZdYrTX6tSHGZO73ezuAw2uxuQ+BbeZfMxsexY35+/Qe4J1inV0g4x3WHuMhIl2jXOKT2Q+VU5fgoSUI52VfZyfkJGkha6jmY+aaYdirq+lQAACH5BAkBAAEALAAAAAAoACgAAAKWBIKpYe0L3YNKToqswUlvznigd4wiR4KhZrKt9Upqip61i9E3vMvxRdHlbEFiEXfk9YARYxOZZD6VQ2pUunBmtRXo1Lf8hMVVcNl8JafV38aM2/Fu5V16Bn63r6xt97j09+MXSFi4BniGFae3hzbH9+hYBzkpuUh5aZmHuanZOZgIuvbGiNeomCnaxxap2upaCZsq+1kAACH5BAkBAAEALAAAAAAoACgAAAKXjI8By5zf4kOxTVrXNVlv1X0d8IGZGKLnNpYtm8Lr9cqVeuOSvfOW79D9aDHizNhDJidFZhNydEahOaDH6nomtJjp1tutKoNWkvA6JqfRVLHU/QUfau9l2x7G54d1fl995xcIGAdXqMfBNadoYrhH+Mg2KBlpVpbluCiXmMnZ2Sh4GBqJ+ckIOqqJ6LmKSllZmsoq6wpQAAAh+QQJAQABACwAAAAAKAAoAAAClYx/oLvoxuJDkU1a1YUZbJ59nSd2ZXhWqbRa2/gF8Gu2DY3iqs7yrq+xBYEkYvFSM8aSSObE+ZgRl1BHFZNr7pRCavZ5BW2142hY3AN/zWtsmf12p9XxxFl2lpLn1rseztfXZjdIWIf2s5dItwjYKBgo9yg5pHgzJXTEeGlZuenpyPmpGQoKOWkYmSpaSnqKileI2FAAACH5BAkBAAEALAAAAAAoACgAAAKVjB+gu+jG4kORTVrVhRlsnn2dJ3ZleFaptFrb+CXmO9OozeL5VfP99HvAWhpiUdcwkpBH3825AwYdU8xTqlLGhtCosArKMpvfa1mMRae9VvWZfeB2XfPkeLmm18lUcBj+p5dnN8jXZ3YIGEhYuOUn45aoCDkp16hl5IjYJvjWKcnoGQpqyPlpOhr3aElaqrq56Bq7VAAAOw==");
	height: 100%;
	filter: alpha(opacity=25); 
	opacity: 0.25;
}
.ui-progressbar-indeterminate .ui-progressbar-value {
	background-image: none;
}
.ui-selectable {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-selectable-helper {
	position: absolute;
	z-index: 100;
	border: 1px dotted black;
}
.ui-selectmenu-menu {
	padding: 0;
	margin: 0;
	position: absolute;
	top: 0;
	left: 0;
	display: none;
}
.ui-selectmenu-menu .ui-menu {
	overflow: auto;
	overflow-x: hidden;
	padding-bottom: 1px;
}
.ui-selectmenu-menu .ui-menu .ui-selectmenu-optgroup {
	font-size: 1em;
	font-weight: bold;
	line-height: 1.5;
	padding: 2px 0.4em;
	margin: 0.5em 0 0 0;
	height: auto;
	border: 0;
}
.ui-selectmenu-open {
	display: block;
}
.ui-selectmenu-text {
	display: block;
	margin-right: 20px;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-selectmenu-button.ui-button {
	text-align: left;
	white-space: nowrap;
	width: 14em;
}
.ui-selectmenu-icon.ui-icon {
	float: right;
	margin-top: 0;
}
.ui-slider {
	position: relative;
	text-align: left;
}
.ui-slider .ui-slider-handle {
	position: absolute;
	z-index: 2;
	width: 1.2em;
	height: 1.2em;
	cursor: default;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-slider .ui-slider-range {
	position: absolute;
	z-index: 1;
	font-size: .7em;
	display: block;
	border: 0;
	background-position: 0 0;
}


.ui-slider.ui-state-disabled .ui-slider-handle,
.ui-slider.ui-state-disabled .ui-slider-range {
	filter: inherit;
}

.ui-slider-horizontal {
	height: .8em;
}
.ui-slider-horizontal .ui-slider-handle {
	top: -.3em;
	margin-left: -.6em;
}
.ui-slider-horizontal .ui-slider-range {
	top: 0;
	height: 100%;
}
.ui-slider-horizontal .ui-slider-range-min {
	left: 0;
}
.ui-slider-horizontal .ui-slider-range-max {
	right: 0;
}

.ui-slider-vertical {
	width: .8em;
	height: 100px;
}
.ui-slider-vertical .ui-slider-handle {
	left: -.3em;
	margin-left: 0;
	margin-bottom: -.6em;
}
.ui-slider-vertical .ui-slider-range {
	left: 0;
	width: 100%;
}
.ui-slider-vertical .ui-slider-range-min {
	bottom: 0;
}
.ui-slider-vertical .ui-slider-range-max {
	top: 0;
}
.ui-sortable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-spinner {
	position: relative;
	display: inline-block;
	overflow: hidden;
	padding: 0;
	vertical-align: middle;
}
.ui-spinner-input {
	border: none;
	background: none;
	color: inherit;
	padding: .222em 0;
	margin: .2em 0;
	vertical-align: middle;
	margin-left: .4em;
	margin-right: 2em;
}
.ui-spinner-button {
	width: 1.6em;
	height: 50%;
	font-size: .5em;
	padding: 0;
	margin: 0;
	text-align: center;
	position: absolute;
	cursor: default;
	display: block;
	overflow: hidden;
	right: 0;
}

.ui-spinner a.ui-spinner-button {
	border-top-style: none;
	border-bottom-style: none;
	border-right-style: none;
}
.ui-spinner-up {
	top: 0;
}
.ui-spinner-down {
	bottom: 0;
}
.ui-tabs {
	position: relative;
	padding: .2em;
}
.ui-tabs .ui-tabs-nav {
	margin: 0;
	padding: .2em .2em 0;
}
.ui-tabs .ui-tabs-nav li {
	list-style: none;
	float: left;
	position: relative;
	top: 0;
	margin: 1px .2em 0 0;
	border-bottom-width: 0;
	padding: 0;
	white-space: nowrap;
}
.ui-tabs .ui-tabs-nav .ui-tabs-anchor {
	float: left;
	padding: .5em 1em;
	text-decoration: none;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active {
	margin-bottom: -1px;
	padding-bottom: 1px;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-state-disabled .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-tabs-loading .ui-tabs-anchor {
	cursor: text;
}
.ui-tabs-collapsible .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor {
	cursor: pointer;
}
.ui-tabs .ui-tabs-panel {
	display: block;
	border-width: 0;
	padding: 1em 1.4em;
	background: none;
}
.ui-tooltip {
	padding: 8px;
	position: absolute;
	z-index: 9999;
	max-width: 300px;
}
body .ui-tooltip {
	border-width: 2px;
}

.ui-widget {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget .ui-widget {
	font-size: 1em;
}
.ui-widget input,
.ui-widget select,
.ui-widget textarea,
.ui-widget button {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget.ui-widget-content {
	border: 1px solid #c5c5c5;
}
.ui-widget-content {
	border: 1px solid #dddddd;
	background: #ffffff;
	color: #333333;
}
.ui-widget-content a {
	color: #333333;
}
.ui-widget-header {
	border: 1px solid #dddddd;
	background: #e9e9e9;
	color: #333333;
	font-weight: bold;
}
.ui-widget-header a {
	color: #333333;
}


.ui-state-default,
.ui-widget-content .ui-state-default,
.ui-widget-header .ui-state-default,
.ui-button,


html .ui-button.ui-state-disabled:hover,
html .ui-button.ui-state-disabled:active {
	border: 1px solid #c5c5c5;
	background: #f6f6f6;
	font-weight: normal;
	color: #454545;
}
.ui-state-default a,
.ui-state-default a:link,
.ui-state-default a:visited,
a.ui-button,
a:link.ui-button,
a:visited.ui-button,
.ui-button {
	color: #454545;
	text-decoration: none;
}
.ui-state-hover,
.ui-widget-content .ui-state-hover,
.ui-widget-header .ui-state-hover,
.ui-state-focus,
.ui-widget-content .ui-state-focus,
.ui-widget-header .ui-state-focus,
.ui-button:hover,
.ui-button:focus {
	border: 1px solid #cccccc;
	background: #ededed;
	font-weight: normal;
	color: #2b2b2b;
}
.ui-state-hover a,
.ui-state-hover a:hover,
.ui-state-hover a:link,
.ui-state-hover a:visited,
.ui-state-focus a,
.ui-state-focus a:hover,
.ui-state-focus a:link,
.ui-state-focus a:visited,
a.ui-button:hover,
a.ui-button:focus {
	color: #2b2b2b;
	text-decoration: none;
}

.ui-visual-focus {
	box-shadow: 0 0 3px 1px rgb(94, 158, 214);
}
.ui-state-active,
.ui-widget-content .ui-state-active,
.ui-widget-header .ui-state-active,
a.ui-button:active,
.ui-button:active,
.ui-button.ui-state-active:hover {
	border: 1px solid #003eff;
	background: #007fff;
	font-weight: normal;
	color: #ffffff;
}
.ui-icon-background,
.ui-state-active .ui-icon-background {
	border: #003eff;
	background-color: #ffffff;
}
.ui-state-active a,
.ui-state-active a:link,
.ui-state-active a:visited {
	color: #ffffff;
	text-decoration: none;
}


.ui-state-highlight,
.ui-widget-content .ui-state-highlight,
.ui-widget-header .ui-state-highlight {
	border: 1px solid #dad55e;
	background: #fffa90;
	color: #777620;
}
.ui-state-checked {
	border: 1px solid #dad55e;
	background: #fffa90;
}
.ui-state-highlight a,
.ui-widget-content .ui-state-highlight a,
.ui-widget-header .ui-state-highlight a {
	color: #777620;
}
.ui-state-error,
.ui-widget-content .ui-state-error,
.ui-widget-header .ui-state-error {
	border: 1px solid #f1a899;
	background: #fddfdf;
	color: #5f3f3f;
}
.ui-state-error a,
.ui-widget-content .ui-state-error a,
.ui-widget-header .ui-state-error a {
	color: #5f3f3f;
}
.ui-state-error-text,
.ui-widget-content .ui-state-error-text,
.ui-widget-header .ui-state-error-text {
	color: #5f3f3f;
}
.ui-priority-primary,
.ui-widget-content .ui-priority-primary,
.ui-widget-header .ui-priority-primary {
	font-weight: bold;
}
.ui-priority-secondary,
.ui-widget-content .ui-priority-secondary,
.ui-widget-header .ui-priority-secondary {
	opacity: .7;
	filter:Alpha(Opacity=70); 
	font-weight: normal;
}
.ui-state-disabled,
.ui-widget-content .ui-state-disabled,
.ui-widget-header .ui-state-disabled {
	opacity: .35;
	filter:Alpha(Opacity=35); 
	background-image: none;
}
.ui-state-disabled .ui-icon {
	filter:Alpha(Opacity=35); 
}




.ui-icon {
	width: 16px;
	height: 16px;
}
.ui-icon,
.ui-widget-content .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-widget-header .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-state-hover .ui-icon,
.ui-state-focus .ui-icon,
.ui-button:hover .ui-icon,
.ui-button:focus .ui-icon {
	background-image: url("images/ui-icons_555555_256x240.png");
}
.ui-state-active .ui-icon,
.ui-button:active .ui-icon {
	background-image: url("images/ui-icons_ffffff_256x240.png");
}
.ui-state-highlight .ui-icon,
.ui-button .ui-state-highlight.ui-icon {
	background-image: url("images/ui-icons_777620_256x240.png");
}
.ui-state-error .ui-icon,
.ui-state-error-text .ui-icon {
	background-image: url("images/ui-icons_cc0000_256x240.png");
}
.ui-button .ui-icon {
	background-image: url("images/ui-icons_777777_256x240.png");
}


.ui-icon-blank { background-position: 16px 16px; }
.ui-icon-caret-1-n { background-position: 0 0; }
.ui-icon-caret-1-ne { background-position: -16px 0; }
.ui-icon-caret-1-e { background-position: -32px 0; }
.ui-icon-caret-1-se { background-position: -48px 0; }
.ui-icon-caret-1-s { background-position: -65px 0; }
.ui-icon-caret-1-sw { background-position: -80px 0; }
.ui-icon-caret-1-w { background-position: -96px 0; }
.ui-icon-caret-1-nw { background-position: -112px 0; }
.ui-icon-caret-2-n-s { background-position: -128px 0; }
.ui-icon-caret-2-e-w { background-position: -144px 0; }
.ui-icon-triangle-1-n { background-position: 0 -16px; }
.ui-icon-triangle-1-ne { background-position: -16px -16px; }
.ui-icon-triangle-1-e { background-position: -32px -16px; }
.ui-icon-triangle-1-se { background-position: -48px -16px; }
.ui-icon-triangle-1-s { background-position: -65px -16px; }
.ui-icon-triangle-1-sw { background-position: -80px -16px; }
.ui-icon-triangle-1-w { background-position: -96px -16px; }
.ui-icon-triangle-1-nw { background-position: -112px -16px; }
.ui-icon-triangle-2-n-s { background-position: -128px -16px; }
.ui-icon-triangle-2-e-w { background-position: -144px -16px; }
.ui-icon-arrow-1-n { background-position: 0 -32px; }
.ui-icon-arrow-1-ne { background-position: -16px -32px; }
.ui-icon-arrow-1-e { background-position: -32px -32px; }
.ui-icon-arrow-1-se { background-position: -48px -32px; }
.ui-icon-arrow-1-s { background-position: -65px -32px; }
.ui-icon-arrow-1-sw { background-position: -80px -32px; }
.ui-icon-arrow-1-w { background-position: -96px -32px; }
.ui-icon-arrow-1-nw { background-position: -112px -32px; }
.ui-icon-arrow-2-n-s { background-position: -128px -32px; }
.ui-icon-arrow-2-ne-sw { background-position: -144px -32px; }
.ui-icon-arrow-2-e-w { background-position: -160px -32px; }
.ui-icon-arrow-2-se-nw { background-position: -176px -32px; }
.ui-icon-arrowstop-1-n { background-position: -192px -32px; }
.ui-icon-arrowstop-1-e { background-position: -208px -32px; }
.ui-icon-arrowstop-1-s { background-position: -224px -32px; }
.ui-icon-arrowstop-1-w { background-position: -240px -32px; }
.ui-icon-arrowthick-1-n { background-position: 1px -48px; }
.ui-icon-arrowthick-1-ne { background-position: -16px -48px; }
.ui-icon-arrowthick-1-e { background-position: -32px -48px; }
.ui-icon-arrowthick-1-se { background-position: -48px -48px; }
.ui-icon-arrowthick-1-s { background-position: -64px -48px; }
.ui-icon-arrowthick-1-sw { background-position: -80px -48px; }
.ui-icon-arrowthick-1-w { background-position: -96px -48px; }
.ui-icon-arrowthick-1-nw { background-position: -112px -48px; }
.ui-icon-arrowthick-2-n-s { background-position: -128px -48px; }
.ui-icon-arrowthick-2-ne-sw { background-position: -144px -48px; }
.ui-icon-arrowthick-2-e-w { background-position: -160px -48px; }
.ui-icon-arrowthick-2-se-nw { background-position: -176px -48px; }
.ui-icon-arrowthickstop-1-n { background-position: -192px -48px; }
.ui-icon-arrowthickstop-1-e { background-position: -208px -48px; }
.ui-icon-arrowthickstop-1-s { background-position: -224px -48px; }
.ui-icon-arrowthickstop-1-w { background-position: -240px -48px; }
.ui-icon-arrowreturnthick-1-w { background-position: 0 -64px; }
.ui-icon-arrowreturnthick-1-n { background-position: -16px -64px; }
.ui-icon-arrowreturnthick-1-e { background-position: -32px -64px; }
.ui-icon-arrowreturnthick-1-s { background-position: -48px -64px; }
.ui-icon-arrowreturn-1-w { background-position: -64px -64px; }
.ui-icon-arrowreturn-1-n { background-position: -80px -64px; }
.ui-icon-arrowreturn-1-e { background-position: -96px -64px; }
.ui-icon-arrowreturn-1-s { background-position: -112px -64px; }
.ui-icon-arrowrefresh-1-w { background-position: -128px -64px; }
.ui-icon-arrowrefresh-1-n { background-position: -144px -64px; }
.ui-icon-arrowrefresh-1-e { background-position: -160px -64px; }
.ui-icon-arrowrefresh-1-s { background-position: -176px -64px; }
.ui-icon-arrow-4 { background-position: 0 -80px; }
.ui-icon-arrow-4-diag { background-position: -16px -80px; }
.ui-icon-extlink { background-position: -32px -80px; }
.ui-icon-newwin { background-position: -48px -80px; }
.ui-icon-refresh { background-position: -64px -80px; }
.ui-icon-shuffle { background-position: -80px -80px; }
.ui-icon-transfer-e-w { background-position: -96px -80px; }
.ui-icon-transferthick-e-w { background-position: -112px -80px; }
.ui-icon-folder-collapsed { background-position: 0 -96px; }
.ui-icon-folder-open { background-position: -16px -96px; }
.ui-icon-document { background-position: -32px -96px; }
.ui-icon-document-b { background-position: -48px -96px; }
.ui-icon-note { background-position: -64px -96px; }
.ui-icon-mail-closed { background-position: -80px -96px; }
.ui-icon-mail-open { background-position: -96px -96px; }
.ui-icon-suitcase { background-position: -112px -96px; }
.ui-icon-comment { background-position: -128px -96px; }
.ui-icon-person { background-position: -144px -96px; }
.ui-icon-print { background-position: -160px -96px; }
.ui-icon-trash { background-position: -176px -96px; }
.ui-icon-locked { background-position: -192px -96px; }
.ui-icon-unlocked { background-position: -208px -96px; }
.ui-icon-bookmark { background-position: -224px -96px; }
.ui-icon-tag { background-position: -240px -96px; }
.ui-icon-home { background-position: 0 -112px; }
.ui-icon-flag { background-position: -16px -112px; }
.ui-icon-calendar { background-position: -32px -112px; }
.ui-icon-cart { background-position: -48px -112px; }
.ui-icon-pencil { background-position: -64px -112px; }
.ui-icon-clock { background-position: -80px -112px; }
.ui-icon-disk { background-position: -96px -112px; }
.ui-icon-calculator { background-position: -112px -112px; }
.ui-icon-zoomin { background-position: -128px -112px; }
.ui-icon-zoomout { background-position: -144px -112px; }
.ui-icon-search { background-position: -160px -112px; }
.ui-icon-wrench { background-position: -176px -112px; }
.ui-icon-gear { background-position: -192px -112px; }
.ui-icon-heart { background-position: -208px -112px; }
.ui-icon-star { background-position: -224px -112px; }
.ui-icon-link { background-position: -240px -112px; }
.ui-icon-cancel { background-position: 0 -128px; }
.ui-icon-plus { background-position: -16px -128px; }
.ui-icon-plusthick { background-position: -32px -128px; }
.ui-icon-minus { background-position: -48px -128px; }
.ui-icon-minusthick { background-position: -64px -128px; }
.ui-icon-close { background-position: -80px -128px; }
.ui-icon-closethick { background-position: -96px -128px; }
.ui-icon-key { background-position: -112px -128px; }
.ui-icon-lightbulb { background-position: -128px -128px; }
.ui-icon-scissors { background-position: -144px -128px; }
.ui-icon-clipboard { background-position: -160px -128px; }
.ui-icon-copy { background-position: -176px -128px; }
.ui-icon-contact { background-position: -192px -128px; }
.ui-icon-image { background-position: -208px -128px; }
.ui-icon-video { background-position: -224px -128px; }
.ui-icon-script { background-position: -240px -128px; }
.ui-icon-alert { background-position: 0 -144px; }
.ui-icon-info { background-position: -16px -144px; }
.ui-icon-notice { background-position: -32px -144px; }
.ui-icon-help { background-position: -48px -144px; }
.ui-icon-check { background-position: -64px -144px; }
.ui-icon-bullet { background-position: -80px -144px; }
.ui-icon-radio-on { background-position: -96px -144px; }
.ui-icon-radio-off { background-position: -112px -144px; }
.ui-icon-pin-w { background-position: -128px -144px; }
.ui-icon-pin-s { background-position: -144px -144px; }
.ui-icon-play { background-position: 0 -160px; }
.ui-icon-pause { background-position: -16px -160px; }
.ui-icon-seek-next { background-position: -32px -160px; }
.ui-icon-seek-prev { background-position: -48px -160px; }
.ui-icon-seek-end { background-position: -64px -160px; }
.ui-icon-seek-start { background-position: -80px -160px; }

.ui-icon-seek-first { background-position: -80px -160px; }
.ui-icon-stop { background-position: -96px -160px; }
.ui-icon-eject { background-position: -112px -160px; }
.ui-icon-volume-off { background-position: -128px -160px; }
.ui-icon-volume-on { background-position: -144px -160px; }
.ui-icon-power { background-position: 0 -176px; }
.ui-icon-signal-diag { background-position: -16px -176px; }
.ui-icon-signal { background-position: -32px -176px; }
.ui-icon-battery-0 { background-position: -48px -176px; }
.ui-icon-battery-1 { background-position: -64px -176px; }
.ui-icon-battery-2 { background-position: -80px -176px; }
.ui-icon-battery-3 { background-position: -96px -176px; }
.ui-icon-circle-plus { background-position: 0 -192px; }
.ui-icon-circle-minus { background-position: -16px -192px; }
.ui-icon-circle-close { background-position: -32px -192px; }
.ui-icon-circle-triangle-e { background-position: -48px -192px; }
.ui-icon-circle-triangle-s { background-position: -64px -192px; }
.ui-icon-circle-triangle-w { background-position: -80px -192px; }
.ui-icon-circle-triangle-n { background-position: -96px -192px; }
.ui-icon-circle-arrow-e { background-position: -112px -192px; }
.ui-icon-circle-arrow-s { background-position: -128px -192px; }
.ui-icon-circle-arrow-w { background-position: -144px -192px; }
.ui-icon-circle-arrow-n { background-position: -160px -192px; }
.ui-icon-circle-zoomin { background-position: -176px -192px; }
.ui-icon-circle-zoomout { background-position: -192px -192px; }
.ui-icon-circle-check { background-position: -208px -192px; }
.ui-icon-circlesmall-plus { background-position: 0 -208px; }
.ui-icon-circlesmall-minus { background-position: -16px -208px; }
.ui-icon-circlesmall-close { background-position: -32px -208px; }
.ui-icon-squaresmall-plus { background-position: -48px -208px; }
.ui-icon-squaresmall-minus { background-position: -64px -208px; }
.ui-icon-squaresmall-close { background-position: -80px -208px; }
.ui-icon-grip-dotted-vertical { background-position: 0 -224px; }
.ui-icon-grip-dotted-horizontal { background-position: -16px -224px; }
.ui-icon-grip-solid-vertical { background-position: -32px -224px; }
.ui-icon-grip-solid-horizontal { background-position: -48px -224px; }
.ui-icon-gripsmall-diagonal-se { background-position: -64px -224px; }
.ui-icon-grip-diagonal-se { background-position: -80px -224px; }





.ui-corner-all,
.ui-corner-top,
.ui-corner-left,
.ui-corner-tl {
	border-top-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-top,
.ui-corner-right,
.ui-corner-tr {
	border-top-right-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-left,
.ui-corner-bl {
	border-bottom-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-right,
.ui-corner-br {
	border-bottom-right-radius: 3px;
}


.ui-widget-overlay {
	background: #aaaaaa;
	opacity: .3;
	filter: Alpha(Opacity=30); 
}
.ui-widget-shadow {
	-webkit-box-shadow: 0px 0px 5px #666666;
	box-shadow: 0px 0px 5px #666666;
} 
</style>
<script src='js/jquery-1.12.4.js'></script>
<script src='js/jquery-ui.js'></script>

<script>
$( function() {
	$( '#tabs-collect-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-collect-consensus-heatmap'>
<ul>
<li><a href='#tab-collect-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-collect-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-collect-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-collect-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-collect-consensus-heatmap-5'>k = 6</a></li>
<li><a href='#tab-collect-consensus-heatmap-6'>k = 7</a></li>
<li><a href='#tab-collect-consensus-heatmap-7'>k = 8</a></li>
<li><a href='#tab-collect-consensus-heatmap-8'>k = 9</a></li>
<li><a href='#tab-collect-consensus-heatmap-9'>k = 10</a></li>
</ul>
<div id='tab-collect-consensus-heatmap-1'>
<pre><code class="r">collect_plots(res_list, k = 2, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-1-1.png" alt="plot of chunk tab-collect-consensus-heatmap-1"/></p>

</div>
<div id='tab-collect-consensus-heatmap-2'>
<pre><code class="r">collect_plots(res_list, k = 3, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-2-1.png" alt="plot of chunk tab-collect-consensus-heatmap-2"/></p>

</div>
<div id='tab-collect-consensus-heatmap-3'>
<pre><code class="r">collect_plots(res_list, k = 4, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-3-1.png" alt="plot of chunk tab-collect-consensus-heatmap-3"/></p>

</div>
<div id='tab-collect-consensus-heatmap-4'>
<pre><code class="r">collect_plots(res_list, k = 5, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-4-1.png" alt="plot of chunk tab-collect-consensus-heatmap-4"/></p>

</div>
<div id='tab-collect-consensus-heatmap-5'>
<pre><code class="r">collect_plots(res_list, k = 6, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-5-1.png" alt="plot of chunk tab-collect-consensus-heatmap-5"/></p>

</div>
<div id='tab-collect-consensus-heatmap-6'>
<pre><code class="r">collect_plots(res_list, k = 7, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-6-1.png" alt="plot of chunk tab-collect-consensus-heatmap-6"/></p>

</div>
<div id='tab-collect-consensus-heatmap-7'>
<pre><code class="r">collect_plots(res_list, k = 8, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-7-1.png" alt="plot of chunk tab-collect-consensus-heatmap-7"/></p>

</div>
<div id='tab-collect-consensus-heatmap-8'>
<pre><code class="r">collect_plots(res_list, k = 9, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-8-1.png" alt="plot of chunk tab-collect-consensus-heatmap-8"/></p>

</div>
<div id='tab-collect-consensus-heatmap-9'>
<pre><code class="r">collect_plots(res_list, k = 10, fun = consensus_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-consensus-heatmap-9-1.png" alt="plot of chunk tab-collect-consensus-heatmap-9"/></p>

</div>
</div>



### Membership heatmap

Membership heatmaps for all methods. ([What is a membership heatmap?](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_12))


<script>
$( function() {
	$( '#tabs-collect-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-collect-membership-heatmap'>
<ul>
<li><a href='#tab-collect-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-collect-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-collect-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-collect-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-collect-membership-heatmap-5'>k = 6</a></li>
<li><a href='#tab-collect-membership-heatmap-6'>k = 7</a></li>
<li><a href='#tab-collect-membership-heatmap-7'>k = 8</a></li>
<li><a href='#tab-collect-membership-heatmap-8'>k = 9</a></li>
<li><a href='#tab-collect-membership-heatmap-9'>k = 10</a></li>
</ul>
<div id='tab-collect-membership-heatmap-1'>
<pre><code class="r">collect_plots(res_list, k = 2, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-1-1.png" alt="plot of chunk tab-collect-membership-heatmap-1"/></p>

</div>
<div id='tab-collect-membership-heatmap-2'>
<pre><code class="r">collect_plots(res_list, k = 3, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-2-1.png" alt="plot of chunk tab-collect-membership-heatmap-2"/></p>

</div>
<div id='tab-collect-membership-heatmap-3'>
<pre><code class="r">collect_plots(res_list, k = 4, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-3-1.png" alt="plot of chunk tab-collect-membership-heatmap-3"/></p>

</div>
<div id='tab-collect-membership-heatmap-4'>
<pre><code class="r">collect_plots(res_list, k = 5, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-4-1.png" alt="plot of chunk tab-collect-membership-heatmap-4"/></p>

</div>
<div id='tab-collect-membership-heatmap-5'>
<pre><code class="r">collect_plots(res_list, k = 6, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-5-1.png" alt="plot of chunk tab-collect-membership-heatmap-5"/></p>

</div>
<div id='tab-collect-membership-heatmap-6'>
<pre><code class="r">collect_plots(res_list, k = 7, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-6-1.png" alt="plot of chunk tab-collect-membership-heatmap-6"/></p>

</div>
<div id='tab-collect-membership-heatmap-7'>
<pre><code class="r">collect_plots(res_list, k = 8, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-7-1.png" alt="plot of chunk tab-collect-membership-heatmap-7"/></p>

</div>
<div id='tab-collect-membership-heatmap-8'>
<pre><code class="r">collect_plots(res_list, k = 9, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-8-1.png" alt="plot of chunk tab-collect-membership-heatmap-8"/></p>

</div>
<div id='tab-collect-membership-heatmap-9'>
<pre><code class="r">collect_plots(res_list, k = 10, fun = membership_heatmap, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-membership-heatmap-9-1.png" alt="plot of chunk tab-collect-membership-heatmap-9"/></p>

</div>
</div>



### Signature heatmap

Signature heatmaps for all methods. ([What is a signature heatmap?](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_21))


Note in following heatmaps, rows are scaled.



<script>
$( function() {
	$( '#tabs-collect-get-signatures' ).tabs();
} );
</script>
<div id='tabs-collect-get-signatures'>
<ul>
<li><a href='#tab-collect-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-collect-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-collect-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-collect-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-collect-get-signatures-5'>k = 6</a></li>
<li><a href='#tab-collect-get-signatures-6'>k = 7</a></li>
<li><a href='#tab-collect-get-signatures-7'>k = 8</a></li>
<li><a href='#tab-collect-get-signatures-8'>k = 9</a></li>
<li><a href='#tab-collect-get-signatures-9'>k = 10</a></li>
</ul>
<div id='tab-collect-get-signatures-1'>
<pre><code class="r">collect_plots(res_list, k = 2, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-1-1.png" alt="plot of chunk tab-collect-get-signatures-1"/></p>

</div>
<div id='tab-collect-get-signatures-2'>
<pre><code class="r">collect_plots(res_list, k = 3, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-2-1.png" alt="plot of chunk tab-collect-get-signatures-2"/></p>

</div>
<div id='tab-collect-get-signatures-3'>
<pre><code class="r">collect_plots(res_list, k = 4, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-3-1.png" alt="plot of chunk tab-collect-get-signatures-3"/></p>

</div>
<div id='tab-collect-get-signatures-4'>
<pre><code class="r">collect_plots(res_list, k = 5, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-4-1.png" alt="plot of chunk tab-collect-get-signatures-4"/></p>

</div>
<div id='tab-collect-get-signatures-5'>
<pre><code class="r">collect_plots(res_list, k = 6, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-5-1.png" alt="plot of chunk tab-collect-get-signatures-5"/></p>

</div>
<div id='tab-collect-get-signatures-6'>
<pre><code class="r">collect_plots(res_list, k = 7, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-6-1.png" alt="plot of chunk tab-collect-get-signatures-6"/></p>

</div>
<div id='tab-collect-get-signatures-7'>
<pre><code class="r">collect_plots(res_list, k = 8, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-7-1.png" alt="plot of chunk tab-collect-get-signatures-7"/></p>

</div>
<div id='tab-collect-get-signatures-8'>
<pre><code class="r">collect_plots(res_list, k = 9, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-8-1.png" alt="plot of chunk tab-collect-get-signatures-8"/></p>

</div>
<div id='tab-collect-get-signatures-9'>
<pre><code class="r">collect_plots(res_list, k = 10, fun = get_signatures, mc.cores = 1)
</code></pre>

<p><img src="figure_cola/tab-collect-get-signatures-9-1.png" alt="plot of chunk tab-collect-get-signatures-9"/></p>

</div>
</div>



### Statistics table

The statistics used for measuring the stability of consensus partitioning.
([How are they
defined?](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_13))


<script>
$( function() {
	$( '#tabs-get-stats-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-get-stats-from-consensus-partition-list'>
<ul>
<li><a href='#tab-get-stats-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-5'>k = 6</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-6'>k = 7</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-7'>k = 8</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-8'>k = 9</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-9'>k = 10</a></li>
</ul>
<div id='tab-get-stats-from-consensus-partition-list-1'>
<pre><code class="r">get_stats(res_list, k = 2)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 2 0.679           0.858       0.937          0.463 0.537   0.537
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-2'>
<pre><code class="r">get_stats(res_list, k = 3)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 3 0.762           0.851       0.931          0.399 0.732   0.534
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-3'>
<pre><code class="r">get_stats(res_list, k = 4)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 4 0.687           0.732       0.834         0.0961 0.898    0.72
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-4'>
<pre><code class="r">get_stats(res_list, k = 5)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 5 0.768           0.783       0.893         0.0656 0.932   0.768
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-5'>
<pre><code class="r">get_stats(res_list, k = 6)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 6  0.76           0.707       0.839         0.0399 0.944   0.777
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-6'>
<pre><code class="r">get_stats(res_list, k = 7)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased Rand Jaccard
#&gt; ATC:skmeans 7 0.759           0.655       0.788         0.0271 0.94   0.737
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-7'>
<pre><code class="r">get_stats(res_list, k = 8)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 8  0.77           0.661       0.784         0.0209 0.963   0.811
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-8'>
<pre><code class="r">get_stats(res_list, k = 9)
</code></pre>

<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 9 0.766           0.608       0.775         0.0193 0.968   0.823
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-9'>
<pre><code class="r">get_stats(res_list, k = 10)
</code></pre>

<pre><code>#&gt;              k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; ATC:skmeans 10 0.778           0.565       0.746         0.0185 0.972   0.828
</code></pre>

</div>
</div>

Following heatmap plots the partition for each combination of methods and the
lightness correspond to the silhouette scores for samples in each method. On
top the consensus subgroup is inferred from all methods by taking the mean
silhouette scores as weight.


<script>
$( function() {
	$( '#tabs-collect-stats-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-collect-stats-from-consensus-partition-list'>
<ul>
<li><a href='#tab-collect-stats-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-5'>k = 6</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-6'>k = 7</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-7'>k = 8</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-8'>k = 9</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-9'>k = 10</a></li>
</ul>
<div id='tab-collect-stats-from-consensus-partition-list-1'>
<pre><code class="r">collect_stats(res_list, k = 2)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-1-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-1"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-2'>
<pre><code class="r">collect_stats(res_list, k = 3)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-2-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-2"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-3'>
<pre><code class="r">collect_stats(res_list, k = 4)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-3-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-3"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-4'>
<pre><code class="r">collect_stats(res_list, k = 5)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-4-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-4"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-5'>
<pre><code class="r">collect_stats(res_list, k = 6)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-5-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-5"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-6'>
<pre><code class="r">collect_stats(res_list, k = 7)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-6-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-6"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-7'>
<pre><code class="r">collect_stats(res_list, k = 8)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-7-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-7"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-8'>
<pre><code class="r">collect_stats(res_list, k = 9)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-8-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-8"/></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-9'>
<pre><code class="r">collect_stats(res_list, k = 10)
</code></pre>

<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-9-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-9"/></p>

</div>
</div>

### Partition from all methods



Collect partitions from all methods:


<script>
$( function() {
	$( '#tabs-collect-classes-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-collect-classes-from-consensus-partition-list'>
<ul>
<li><a href='#tab-collect-classes-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-5'>k = 6</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-6'>k = 7</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-7'>k = 8</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-8'>k = 9</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-9'>k = 10</a></li>
</ul>
<div id='tab-collect-classes-from-consensus-partition-list-1'>
<pre><code class="r">collect_classes(res_list, k = 2)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-1-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-1"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-2'>
<pre><code class="r">collect_classes(res_list, k = 3)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-2-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-2"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-3'>
<pre><code class="r">collect_classes(res_list, k = 4)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-3-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-3"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-4'>
<pre><code class="r">collect_classes(res_list, k = 5)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-4-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-4"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-5'>
<pre><code class="r">collect_classes(res_list, k = 6)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-5-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-5"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-6'>
<pre><code class="r">collect_classes(res_list, k = 7)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-6-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-6"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-7'>
<pre><code class="r">collect_classes(res_list, k = 8)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-7-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-7"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-8'>
<pre><code class="r">collect_classes(res_list, k = 9)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-8-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-8"/></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-9'>
<pre><code class="r">collect_classes(res_list, k = 10)
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-9-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-9"/></p>

</div>
</div>



### Top rows overlap


Heatmaps of the top rows:



<script>
$( function() {
	$( '#tabs-top-rows-heatmap' ).tabs();
} );
</script>
<div id='tabs-top-rows-heatmap'>
<ul>
<li><a href='#tab-top-rows-heatmap-1'>top_n = 13</a></li>
</ul>
<div id='tab-top-rows-heatmap-1'>
<pre><code class="r">top_rows_heatmap(res_list, top_n = 13)
</code></pre>

<p><img src="figure_cola/tab-top-rows-heatmap-1-1.png" alt="plot of chunk tab-top-rows-heatmap-1"/></p>

</div>
</div>



 
## Results for each method


---------------------------------------------------



### ATC:skmeans






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_list["ATC", "skmeans"]
# you can also extract it by
# res = res_list["ATC:skmeans"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6, 7, 8, 9, 10.
#>   On a matrix with 13 rows and 2737 columns.
#>   Top rows (13) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 450 partitions by row resampling.
#>   Best k for subgroups seems to be 2.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_signatures"     
#>  [7] "consensus_heatmap"       "dimension_reduction"     "functional_enrichment"  
#> [10] "get_anno_col"            "get_anno"                "get_classes"            
#> [13] "get_consensus"           "get_matrix"              "get_membership"         
#> [16] "get_param"               "get_signatures"          "get_stats"              
#> [19] "is_best_k"               "is_stable_k"             "membership_heatmap"     
#> [22] "ncol"                    "nrow"                    "plot_ecdf"              
#> [25] "predict_classes"         "rownames"                "select_partition_number"
#> [28] "show"                    "suggest_best_k"          "test_to_known_factors"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk ATC-skmeans-collect-plots](figure_cola/ATC-skmeans-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk ATC-skmeans-select-partition-number](figure_cola/ATC-skmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>     k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2   2 0.679           0.858       0.937         0.4626 0.537   0.537
#> 3   3 0.762           0.851       0.931         0.3990 0.732   0.534
#> 4   4 0.687           0.732       0.834         0.0961 0.898   0.720
#> 5   5 0.768           0.783       0.893         0.0656 0.932   0.768
#> 6   6 0.760           0.707       0.839         0.0399 0.944   0.777
#> 7   7 0.759           0.655       0.788         0.0271 0.940   0.737
#> 8   8 0.770           0.661       0.784         0.0209 0.963   0.811
#> 9   9 0.766           0.608       0.775         0.0193 0.968   0.823
#> 10 10 0.778           0.565       0.746         0.0185 0.972   0.828
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 2
```


Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-classes'>
<ul>
<li><a href='#tab-ATC-skmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-9'>k = 10</a></li>
</ul>

<div id='tab-ATC-skmeans-get-classes-1'>
<p><a id='tab-ATC-skmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2
#&gt; SP1003       2   0.000    0.90780 0.00 1.00
#&gt; SP10084      2   0.000    0.90780 0.00 1.00
#&gt; SP1009       1   0.904    0.55054 0.68 0.32
#&gt; SP10150      2   0.327    0.88640 0.06 0.94
#&gt; SP101515     2   0.855    0.62118 0.28 0.72
#&gt; SP101519     2   0.000    0.90780 0.00 1.00
#&gt; SP101521     2   0.000    0.90780 0.00 1.00
#&gt; SP101523     1   0.827    0.66400 0.74 0.26
#&gt; SP101526     2   0.000    0.90780 0.00 1.00
#&gt; SP101528     2   0.529    0.83826 0.12 0.88
#&gt; SP101532     1   0.680    0.78093 0.82 0.18
#&gt; SP101536     2   0.242    0.89595 0.04 0.96
#&gt; SP101540     2   0.000    0.90780 0.00 1.00
#&gt; SP101544     2   0.000    0.90780 0.00 1.00
#&gt; SP101548     2   0.999    0.09622 0.48 0.52
#&gt; SP101552     2   0.000    0.90780 0.00 1.00
#&gt; SP101558     1   0.760    0.72593 0.78 0.22
#&gt; SP101564     2   0.000    0.90780 0.00 1.00
#&gt; SP101572     2   0.141    0.90316 0.02 0.98
#&gt; SP101576     2   0.000    0.90780 0.00 1.00
#&gt; SP101580     2   0.000    0.90780 0.00 1.00
#&gt; SP101584     2   0.000    0.90780 0.00 1.00
#&gt; SP101588     2   0.000    0.90780 0.00 1.00
#&gt; SP101592     2   0.000    0.90780 0.00 1.00
#&gt; SP101596     2   0.000    0.90780 0.00 1.00
#&gt; SP101600     1   0.795    0.69454 0.76 0.24
#&gt; SP101604     1   0.827    0.66257 0.74 0.26
#&gt; SP101610     1   0.990    0.23082 0.56 0.44
#&gt; SP101616     2   0.000    0.90780 0.00 1.00
#&gt; SP101622     2   0.000    0.90780 0.00 1.00
#&gt; SP101628     2   0.000    0.90780 0.00 1.00
#&gt; SP101634     1   0.943    0.45989 0.64 0.36
#&gt; SP101642     2   0.000    0.90780 0.00 1.00
#&gt; SP101648     1   0.999    0.09111 0.52 0.48
#&gt; SP101654     2   0.000    0.90780 0.00 1.00
#&gt; SP101658     1   0.000    0.94439 1.00 0.00
#&gt; SP101662     1   0.469    0.86690 0.90 0.10
#&gt; SP101666     1   0.760    0.72594 0.78 0.22
#&gt; SP101670     2   0.000    0.90780 0.00 1.00
#&gt; SP101674     2   0.000    0.90780 0.00 1.00
#&gt; SP101678     2   0.000    0.90780 0.00 1.00
#&gt; SP101682     2   0.000    0.90780 0.00 1.00
#&gt; SP101686     2   0.000    0.90780 0.00 1.00
#&gt; SP101690     2   0.000    0.90780 0.00 1.00
#&gt; SP101694     2   0.141    0.90316 0.02 0.98
#&gt; SP101700     2   0.141    0.90171 0.02 0.98
#&gt; SP101708     1   0.981    0.29632 0.58 0.42
#&gt; SP101716     2   0.000    0.90780 0.00 1.00
#&gt; SP101724     2   0.000    0.90780 0.00 1.00
#&gt; SP101732     2   0.000    0.90780 0.00 1.00
#&gt; SP101740     2   0.000    0.90780 0.00 1.00
#&gt; SP101795     2   0.000    0.90780 0.00 1.00
#&gt; SP101845     2   0.000    0.90780 0.00 1.00
#&gt; SP101881     1   1.000    0.00885 0.50 0.50
#&gt; SP101891     2   0.000    0.90780 0.00 1.00
#&gt; SP101921     2   0.000    0.90780 0.00 1.00
#&gt; SP101931     2   0.000    0.90780 0.00 1.00
#&gt; SP102015     2   0.000    0.90780 0.00 1.00
#&gt; SP102035     2   0.999    0.07677 0.48 0.52
#&gt; SP102045     2   0.141    0.90316 0.02 0.98
#&gt; SP102055     1   0.722    0.75321 0.80 0.20
#&gt; SP102064     1   0.402    0.88371 0.92 0.08
#&gt; SP102074     2   0.000    0.90780 0.00 1.00
#&gt; SP102084     1   0.971    0.35326 0.60 0.40
#&gt; SP102090     2   0.000    0.90780 0.00 1.00
#&gt; SP102096     2   0.000    0.90780 0.00 1.00
#&gt; SP102103     2   0.000    0.90780 0.00 1.00
#&gt; SP102113     2   0.999    0.06798 0.48 0.52
#&gt; SP102123     2   0.000    0.90780 0.00 1.00
#&gt; SP102133     2   0.000    0.90780 0.00 1.00
#&gt; SP102143     2   0.795    0.69735 0.24 0.76
#&gt; SP102161     2   0.402    0.87577 0.08 0.92
#&gt; SP102168     2   0.680    0.77482 0.18 0.82
#&gt; SP102174     2   0.000    0.90780 0.00 1.00
#&gt; SP102187     1   0.680    0.78120 0.82 0.18
#&gt; SP102485     1   0.000    0.94439 1.00 0.00
#&gt; SP102489     1   0.000    0.94439 1.00 0.00
#&gt; SP102499     1   0.000    0.94439 1.00 0.00
#&gt; SP102507     1   0.402    0.88371 0.92 0.08
#&gt; SP102511     1   0.000    0.94439 1.00 0.00
#&gt; SP102517     1   0.999   -0.06455 0.52 0.48
#&gt; SP102523     1   0.958    0.31009 0.62 0.38
#&gt; SP102529     2   0.999    0.20520 0.48 0.52
#&gt; SP102537     1   0.000    0.94439 1.00 0.00
#&gt; SP102541     2   0.958    0.48163 0.38 0.62
#&gt; SP102547     1   0.000    0.94439 1.00 0.00
#&gt; SP102557     1   0.000    0.94439 1.00 0.00
#&gt; SP102561     1   0.000    0.94439 1.00 0.00
#&gt; SP102567     1   0.990    0.09604 0.56 0.44
#&gt; SP102573     1   0.584    0.80112 0.86 0.14
#&gt; SP102581     1   0.000    0.94439 1.00 0.00
#&gt; SP102591     1   0.000    0.94439 1.00 0.00
#&gt; SP102597     1   0.000    0.94439 1.00 0.00
#&gt; SP102605     1   0.000    0.94439 1.00 0.00
#&gt; SP102611     2   0.242    0.89595 0.04 0.96
#&gt; SP102617     1   0.000    0.94439 1.00 0.00
#&gt; SP102620     1   0.141    0.93033 0.98 0.02
#&gt; SP102622     2   0.958    0.47831 0.38 0.62
#&gt; SP102626     1   0.141    0.93033 0.98 0.02
#&gt; SP102630     1   0.999    0.08349 0.52 0.48
#&gt; SP102633     1   0.000    0.94439 1.00 0.00
#&gt; SP102647     1   0.469    0.86690 0.90 0.10
#&gt; SP102652     2   0.584    0.82076 0.14 0.86
#&gt; SP102690     1   0.000    0.94439 1.00 0.00
#&gt; SP102716     1   0.000    0.94439 1.00 0.00
#&gt; SP102718     1   0.000    0.94439 1.00 0.00
#&gt; SP102733     1   0.000    0.94439 1.00 0.00
#&gt; SP102741     1   0.000    0.94439 1.00 0.00
#&gt; SP102747     1   0.000    0.94439 1.00 0.00
#&gt; SP102755     1   0.000    0.94439 1.00 0.00
#&gt; SP102759     1   0.000    0.94439 1.00 0.00
#&gt; SP102783     1   0.000    0.94439 1.00 0.00
#&gt; SP102804     1   0.000    0.94439 1.00 0.00
#&gt; SP102816     1   0.000    0.94439 1.00 0.00
#&gt; SP102827     2   1.000    0.13899 0.50 0.50
#&gt; SP102839     1   0.000    0.94439 1.00 0.00
#&gt; SP102873     1   0.000    0.94439 1.00 0.00
#&gt; SP102881     1   0.000    0.94439 1.00 0.00
#&gt; SP102897     1   0.000    0.94439 1.00 0.00
#&gt; SP102913     2   0.529    0.85049 0.12 0.88
#&gt; SP102921     1   0.000    0.94439 1.00 0.00
#&gt; SP102929     1   0.000    0.94439 1.00 0.00
#&gt; SP102945     1   0.000    0.94439 1.00 0.00
#&gt; SP102957     1   0.000    0.94439 1.00 0.00
#&gt; SP102965     1   0.000    0.94439 1.00 0.00
#&gt; SP102973     1   0.000    0.94439 1.00 0.00
#&gt; SP102989     1   0.000    0.94439 1.00 0.00
#&gt; SP103005     2   0.855    0.66635 0.28 0.72
#&gt; SP103021     1   0.971    0.23871 0.60 0.40
#&gt; SP103037     1   0.000    0.94439 1.00 0.00
#&gt; SP103045     2   0.795    0.72186 0.24 0.76
#&gt; SP103057     1   0.000    0.94439 1.00 0.00
#&gt; SP103065     1   0.000    0.94439 1.00 0.00
#&gt; SP103080     1   0.795    0.66806 0.76 0.24
#&gt; SP103100     1   0.000    0.94439 1.00 0.00
#&gt; SP103128     1   0.000    0.94439 1.00 0.00
#&gt; SP103140     1   0.000    0.94439 1.00 0.00
#&gt; SP103156     1   0.000    0.94439 1.00 0.00
#&gt; SP103197     1   0.000    0.94439 1.00 0.00
#&gt; SP103213     1   0.000    0.94439 1.00 0.00
#&gt; SP103221     1   0.000    0.94439 1.00 0.00
#&gt; SP103233     1   0.000    0.94439 1.00 0.00
#&gt; SP103245     1   0.000    0.94439 1.00 0.00
#&gt; SP103261     1   0.000    0.94439 1.00 0.00
#&gt; SP103288     1   0.000    0.94439 1.00 0.00
#&gt; SP103300     1   0.000    0.94439 1.00 0.00
#&gt; SP103340     1   0.000    0.94439 1.00 0.00
#&gt; SP103396     1   0.000    0.94439 1.00 0.00
#&gt; SP103408     1   0.000    0.94439 1.00 0.00
#&gt; SP103416     1   0.242    0.91169 0.96 0.04
#&gt; SP103428     1   0.000    0.94439 1.00 0.00
#&gt; SP103436     1   0.000    0.94439 1.00 0.00
#&gt; SP103444     1   0.000    0.94439 1.00 0.00
#&gt; SP103455     2   0.760    0.73804 0.22 0.78
#&gt; SP103467     1   0.000    0.94439 1.00 0.00
#&gt; SP103483     1   0.000    0.94439 1.00 0.00
#&gt; SP103507     1   0.000    0.94439 1.00 0.00
#&gt; SP103515     1   0.000    0.94439 1.00 0.00
#&gt; SP103523     1   0.000    0.94439 1.00 0.00
#&gt; SP103535     2   0.855    0.66598 0.28 0.72
#&gt; SP103547     1   0.000    0.94439 1.00 0.00
#&gt; SP103555     1   0.000    0.94439 1.00 0.00
#&gt; SP103575     2   0.402    0.87577 0.08 0.92
#&gt; SP103595     1   1.000   -0.13730 0.50 0.50
#&gt; SP103603     1   1.000   -0.13859 0.50 0.50
#&gt; SP103619     1   0.000    0.94439 1.00 0.00
#&gt; SP103673     1   0.000    0.94439 1.00 0.00
#&gt; SP103679     1   0.000    0.94439 1.00 0.00
#&gt; SP103685     1   0.000    0.94439 1.00 0.00
#&gt; SP103694     1   0.000    0.94439 1.00 0.00
#&gt; SP103706     1   0.000    0.94439 1.00 0.00
#&gt; SP103715     1   0.000    0.94439 1.00 0.00
#&gt; SP103730     1   0.000    0.94439 1.00 0.00
#&gt; SP103742     1   0.000    0.94439 1.00 0.00
#&gt; SP103826     1   0.000    0.94439 1.00 0.00
#&gt; SP103844     1   0.000    0.94439 1.00 0.00
#&gt; SP103856     1   0.000    0.94439 1.00 0.00
#&gt; SP103866     1   0.000    0.94439 1.00 0.00
#&gt; SP103894     2   0.402    0.87767 0.08 0.92
#&gt; SP104056     1   0.000    0.94439 1.00 0.00
#&gt; SP104330     1   0.000    0.94439 1.00 0.00
#&gt; SP104530     2   0.000    0.90780 0.00 1.00
#&gt; SP10470      2   0.000    0.90780 0.00 1.00
#&gt; SP104984     1   0.584    0.81532 0.86 0.14
#&gt; SP105006     2   0.402    0.87577 0.08 0.92
#&gt; SP105018     2   0.000    0.90780 0.00 1.00
#&gt; SP105086     2   0.529    0.84118 0.12 0.88
#&gt; SP105159     1   0.000    0.94439 1.00 0.00
#&gt; SP105213     1   0.000    0.94439 1.00 0.00
#&gt; SP105253     2   0.000    0.90780 0.00 1.00
#&gt; SP105261     2   0.000    0.90780 0.00 1.00
#&gt; SP105375     2   0.000    0.90780 0.00 1.00
#&gt; SP105425     1   0.827    0.66293 0.74 0.26
#&gt; SP105577     1   0.904    0.54874 0.68 0.32
#&gt; SP10563      2   0.000    0.90780 0.00 1.00
#&gt; SP105673     2   0.999    0.08361 0.48 0.52
#&gt; SP105708     1   0.000    0.94439 1.00 0.00
#&gt; SP105759     1   0.000    0.94439 1.00 0.00
#&gt; SP105807     1   0.000    0.94439 1.00 0.00
#&gt; SP1059       2   0.000    0.90780 0.00 1.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2537 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-1-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-2'>
<p><a id='tab-ATC-skmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3
#&gt; SP1003       2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP10084      2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP1009       3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP10150      2  0.0892     0.9147 0.02 0.98 0.00
#&gt; SP101515     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101519     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101521     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101523     3  0.0892     0.8683 0.02 0.00 0.98
#&gt; SP101526     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101528     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101532     3  0.0892     0.8678 0.02 0.00 0.98
#&gt; SP101536     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101540     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101544     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101548     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101552     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101558     3  0.3340     0.8339 0.12 0.00 0.88
#&gt; SP101564     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101572     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101576     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101580     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101584     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101588     2  0.1529     0.8930 0.00 0.96 0.04
#&gt; SP101592     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101596     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101600     3  0.1529     0.8663 0.04 0.00 0.96
#&gt; SP101604     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101610     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101616     2  0.4291     0.7345 0.00 0.82 0.18
#&gt; SP101622     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101628     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101634     3  0.1529     0.8654 0.04 0.00 0.96
#&gt; SP101642     2  0.0892     0.9097 0.00 0.98 0.02
#&gt; SP101648     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101654     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101658     3  0.5948     0.5244 0.36 0.00 0.64
#&gt; SP101662     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101666     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101670     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101674     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101678     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101682     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101686     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101690     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101694     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101700     3  0.2959     0.8175 0.00 0.10 0.90
#&gt; SP101708     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101716     2  0.6244     0.1654 0.00 0.56 0.44
#&gt; SP101724     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101732     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101740     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101795     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101845     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101881     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP101891     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101921     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP101931     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102015     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102035     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102045     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102055     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102064     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102074     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102084     3  0.0892     0.8683 0.02 0.00 0.98
#&gt; SP102090     2  0.4555     0.6997 0.00 0.80 0.20
#&gt; SP102096     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102103     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102113     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102123     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102133     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102143     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102161     2  0.0892     0.9147 0.02 0.98 0.00
#&gt; SP102168     3  0.7395     0.3871 0.04 0.38 0.58
#&gt; SP102174     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102187     3  0.2537     0.8543 0.08 0.00 0.92
#&gt; SP102485     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102489     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102499     1  0.2537     0.8799 0.92 0.00 0.08
#&gt; SP102507     3  0.2537     0.8603 0.08 0.00 0.92
#&gt; SP102511     1  0.2959     0.8617 0.90 0.00 0.10
#&gt; SP102517     2  0.6045     0.4516 0.38 0.62 0.00
#&gt; SP102523     1  0.6280     0.0653 0.54 0.46 0.00
#&gt; SP102529     2  0.5948     0.4959 0.36 0.64 0.00
#&gt; SP102537     1  0.6309    -0.1001 0.50 0.00 0.50
#&gt; SP102541     2  0.5560     0.6149 0.30 0.70 0.00
#&gt; SP102547     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102557     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102561     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102567     2  0.6302     0.1721 0.48 0.52 0.00
#&gt; SP102573     1  0.3340     0.8281 0.88 0.12 0.00
#&gt; SP102581     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102591     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102597     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102605     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102611     2  0.0000     0.9233 0.00 1.00 0.00
#&gt; SP102617     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102620     3  0.1529     0.8664 0.04 0.00 0.96
#&gt; SP102622     2  0.5216     0.6783 0.26 0.74 0.00
#&gt; SP102626     3  0.1529     0.8664 0.04 0.00 0.96
#&gt; SP102630     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102633     3  0.4796     0.7493 0.22 0.00 0.78
#&gt; SP102647     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102652     3  0.0000     0.8669 0.00 0.00 1.00
#&gt; SP102690     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102716     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102718     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102733     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102741     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102747     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102755     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102759     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102783     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102804     3  0.5560     0.6264 0.30 0.00 0.70
#&gt; SP102816     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102827     2  0.5835     0.5400 0.34 0.66 0.00
#&gt; SP102839     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102873     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102881     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102897     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102913     2  0.3340     0.8373 0.12 0.88 0.00
#&gt; SP102921     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102929     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102945     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102957     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102965     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102973     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP102989     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103005     2  0.5216     0.6792 0.26 0.74 0.00
#&gt; SP103021     1  0.6192     0.2012 0.58 0.42 0.00
#&gt; SP103037     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103045     2  0.4796     0.7272 0.22 0.78 0.00
#&gt; SP103057     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103065     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103080     1  0.9083     0.2043 0.52 0.16 0.32
#&gt; SP103100     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103128     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103140     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103156     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103197     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103213     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103221     3  0.6280     0.2941 0.46 0.00 0.54
#&gt; SP103233     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103245     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103261     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103288     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103300     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103340     1  0.4002     0.7777 0.84 0.00 0.16
#&gt; SP103396     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103408     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103416     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103428     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103436     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103444     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103455     3  0.6407     0.7275 0.08 0.16 0.76
#&gt; SP103467     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103483     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103507     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103515     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103523     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103535     2  0.2537     0.8749 0.08 0.92 0.00
#&gt; SP103547     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103555     1  0.0000     0.9561 1.00 0.00 0.00
#&gt; SP103575     2  0.2066     0.8908 0.06 0.94 0.00
#&gt; SP103595     2  0.6302     0.1721 0.48 0.52 0.00
#&gt; SP103603     2  0.6280     0.2365 0.46 0.54 0.00
#&gt; SP103619     1  0.0000     0.9561 1.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2571 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-2-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-3'>
<p><a id='tab-ATC-skmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4
#&gt; SP1003       2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP10084      2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP1009       4  0.3975     0.6823 0.00 0.00 0.24 0.76
#&gt; SP10150      2  0.3606     0.7611 0.02 0.84 0.00 0.14
#&gt; SP101515     3  0.4948     0.1715 0.00 0.00 0.56 0.44
#&gt; SP101519     2  0.3037     0.7594 0.00 0.88 0.02 0.10
#&gt; SP101521     2  0.6976     0.7339 0.00 0.58 0.18 0.24
#&gt; SP101523     3  0.3335     0.6568 0.02 0.00 0.86 0.12
#&gt; SP101526     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101528     4  0.3975     0.6823 0.00 0.00 0.24 0.76
#&gt; SP101532     3  0.5062     0.4786 0.02 0.00 0.68 0.30
#&gt; SP101536     2  0.3821     0.7671 0.00 0.84 0.04 0.12
#&gt; SP101540     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101544     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101548     4  0.3975     0.6823 0.00 0.00 0.24 0.76
#&gt; SP101552     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101558     3  0.3611     0.6585 0.06 0.00 0.86 0.08
#&gt; SP101564     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101572     2  0.5594     0.7734 0.00 0.72 0.10 0.18
#&gt; SP101576     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101580     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101584     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101588     2  0.5657     0.7719 0.00 0.72 0.12 0.16
#&gt; SP101592     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101596     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101600     3  0.3611     0.6550 0.08 0.00 0.86 0.06
#&gt; SP101604     4  0.4277     0.6311 0.00 0.00 0.28 0.72
#&gt; SP101610     3  0.4624     0.4147 0.00 0.00 0.66 0.34
#&gt; SP101616     3  0.7004     0.0783 0.00 0.20 0.58 0.22
#&gt; SP101622     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101628     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101634     3  0.3335     0.6564 0.02 0.00 0.86 0.12
#&gt; SP101642     2  0.7653     0.5767 0.00 0.46 0.30 0.24
#&gt; SP101648     3  0.2921     0.6463 0.00 0.00 0.86 0.14
#&gt; SP101654     2  0.7583     0.6277 0.00 0.48 0.28 0.24
#&gt; SP101658     3  0.6595     0.5020 0.24 0.04 0.66 0.06
#&gt; SP101662     4  0.4907     0.3437 0.00 0.00 0.42 0.58
#&gt; SP101666     3  0.4522     0.4534 0.00 0.00 0.68 0.32
#&gt; SP101670     2  0.5383     0.7661 0.00 0.74 0.10 0.16
#&gt; SP101674     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101678     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101682     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101686     2  0.5594     0.7728 0.00 0.72 0.10 0.18
#&gt; SP101690     2  0.7394     0.6770 0.00 0.52 0.24 0.24
#&gt; SP101694     2  0.4949     0.7727 0.00 0.76 0.06 0.18
#&gt; SP101700     3  0.3935     0.6462 0.00 0.10 0.84 0.06
#&gt; SP101708     3  0.2921     0.6463 0.00 0.00 0.86 0.14
#&gt; SP101716     3  0.6976     0.0923 0.00 0.18 0.58 0.24
#&gt; SP101724     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101732     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101740     2  0.7135     0.7115 0.00 0.56 0.20 0.24
#&gt; SP101795     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101845     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101881     3  0.2921     0.6463 0.00 0.00 0.86 0.14
#&gt; SP101891     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP101921     2  0.0000     0.7260 0.00 1.00 0.00 0.00
#&gt; SP101931     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102015     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102035     3  0.4522     0.4534 0.00 0.00 0.68 0.32
#&gt; SP102045     2  0.6104     0.7720 0.00 0.68 0.14 0.18
#&gt; SP102055     3  0.4855     0.2924 0.00 0.00 0.60 0.40
#&gt; SP102064     3  0.4522     0.4534 0.00 0.00 0.68 0.32
#&gt; SP102074     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102084     3  0.3037     0.6657 0.02 0.00 0.88 0.10
#&gt; SP102090     3  0.6594     0.2038 0.00 0.14 0.62 0.24
#&gt; SP102096     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102103     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102113     3  0.2921     0.6463 0.00 0.00 0.86 0.14
#&gt; SP102123     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102133     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102143     3  0.4855     0.2681 0.00 0.00 0.60 0.40
#&gt; SP102161     2  0.0707     0.7215 0.02 0.98 0.00 0.00
#&gt; SP102168     3  0.3853     0.5545 0.00 0.02 0.82 0.16
#&gt; SP102174     2  0.6594     0.7638 0.00 0.62 0.14 0.24
#&gt; SP102187     3  0.3611     0.6585 0.06 0.00 0.86 0.08
#&gt; SP102485     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102489     1  0.3172     0.7808 0.84 0.16 0.00 0.00
#&gt; SP102499     1  0.2345     0.8497 0.90 0.00 0.10 0.00
#&gt; SP102507     4  0.5486     0.7088 0.08 0.00 0.20 0.72
#&gt; SP102511     1  0.3037     0.8325 0.88 0.02 0.10 0.00
#&gt; SP102517     2  0.4134     0.5032 0.26 0.74 0.00 0.00
#&gt; SP102523     2  0.4624     0.3993 0.34 0.66 0.00 0.00
#&gt; SP102529     2  0.3801     0.5373 0.22 0.78 0.00 0.00
#&gt; SP102537     3  0.5636     0.4986 0.26 0.00 0.68 0.06
#&gt; SP102541     2  0.3801     0.5373 0.22 0.78 0.00 0.00
#&gt; SP102547     1  0.3172     0.7808 0.84 0.16 0.00 0.00
#&gt; SP102557     1  0.2706     0.8578 0.90 0.02 0.08 0.00
#&gt; SP102561     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102567     2  0.3975     0.5120 0.24 0.76 0.00 0.00
#&gt; SP102573     1  0.4994     0.1609 0.52 0.48 0.00 0.00
#&gt; SP102581     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102591     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102597     1  0.2011     0.8740 0.92 0.00 0.08 0.00
#&gt; SP102605     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102611     2  0.3198     0.7566 0.00 0.88 0.08 0.04
#&gt; SP102617     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102620     4  0.4755     0.7227 0.04 0.00 0.20 0.76
#&gt; SP102622     2  0.2011     0.6906 0.08 0.92 0.00 0.00
#&gt; SP102626     4  0.4755     0.7227 0.04 0.00 0.20 0.76
#&gt; SP102630     4  0.4134     0.6592 0.00 0.00 0.26 0.74
#&gt; SP102633     4  0.5077     0.7421 0.16 0.00 0.08 0.76
#&gt; SP102647     4  0.3975     0.6823 0.00 0.00 0.24 0.76
#&gt; SP102652     4  0.5077     0.6817 0.00 0.08 0.16 0.76
#&gt; SP102690     1  0.0707     0.9362 0.98 0.00 0.00 0.02
#&gt; SP102716     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102718     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102733     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102741     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102747     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102755     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102759     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102783     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102804     4  0.4755     0.7305 0.20 0.00 0.04 0.76
#&gt; SP102816     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102827     2  0.4406     0.4573 0.30 0.70 0.00 0.00
#&gt; SP102839     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102873     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102881     1  0.0707     0.9362 0.98 0.00 0.00 0.02
#&gt; SP102897     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102913     2  0.2011     0.6893 0.08 0.92 0.00 0.00
#&gt; SP102921     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102929     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102945     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102957     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102965     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102973     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP102989     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103005     2  0.2647     0.6518 0.12 0.88 0.00 0.00
#&gt; SP103021     2  0.4406     0.4276 0.30 0.70 0.00 0.00
#&gt; SP103037     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103045     2  0.2345     0.6715 0.10 0.90 0.00 0.00
#&gt; SP103057     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103065     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103080     4  0.7832     0.2704 0.26 0.36 0.00 0.38
#&gt; SP103100     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103128     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103140     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103156     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103197     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103213     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt; SP103221     3  0.4994     0.1912 0.48 0.00 0.52 0.00
#&gt; SP103233     1  0.0000     0.9544 1.00 0.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2595 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-3-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-4'>
<p><a id='tab-ATC-skmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5
#&gt; SP1003       2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP10084      2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP1009       4  0.1043     0.7853 0.00 0.00 0.04 0.96 0.00
#&gt; SP10150      2  0.4287     0.2284 0.00 0.54 0.00 0.00 0.46
#&gt; SP101515     3  0.4262     0.2620 0.00 0.00 0.56 0.44 0.00
#&gt; SP101519     5  0.4060     0.3771 0.00 0.36 0.00 0.00 0.64
#&gt; SP101521     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101523     3  0.0609     0.7715 0.02 0.00 0.98 0.00 0.00
#&gt; SP101526     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101528     4  0.3109     0.6982 0.00 0.00 0.20 0.80 0.00
#&gt; SP101532     3  0.4675     0.5364 0.02 0.00 0.60 0.38 0.00
#&gt; SP101536     2  0.4182     0.4064 0.00 0.60 0.00 0.00 0.40
#&gt; SP101540     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101544     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101548     4  0.0609     0.7944 0.00 0.00 0.02 0.98 0.00
#&gt; SP101552     2  0.0609     0.8695 0.00 0.98 0.00 0.00 0.02
#&gt; SP101558     3  0.0609     0.7715 0.02 0.00 0.98 0.00 0.00
#&gt; SP101564     2  0.0609     0.8693 0.00 0.98 0.00 0.00 0.02
#&gt; SP101572     2  0.2929     0.7559 0.00 0.82 0.00 0.00 0.18
#&gt; SP101576     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101580     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101584     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101588     2  0.4060     0.4593 0.00 0.64 0.00 0.00 0.36
#&gt; SP101592     2  0.0609     0.8693 0.00 0.98 0.00 0.00 0.02
#&gt; SP101596     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101600     3  0.0609     0.7715 0.02 0.00 0.98 0.00 0.00
#&gt; SP101604     4  0.2516     0.7315 0.00 0.00 0.14 0.86 0.00
#&gt; SP101610     3  0.3561     0.6052 0.00 0.00 0.74 0.26 0.00
#&gt; SP101616     2  0.3983     0.5096 0.00 0.66 0.34 0.00 0.00
#&gt; SP101622     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101628     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101634     3  0.1732     0.7610 0.00 0.00 0.92 0.08 0.00
#&gt; SP101642     2  0.1732     0.8318 0.00 0.92 0.08 0.00 0.00
#&gt; SP101648     3  0.0609     0.7690 0.00 0.00 0.98 0.02 0.00
#&gt; SP101654     2  0.0609     0.8681 0.00 0.98 0.02 0.00 0.00
#&gt; SP101658     3  0.6484     0.5681 0.22 0.00 0.60 0.14 0.04
#&gt; SP101662     4  0.4227     0.0702 0.00 0.00 0.42 0.58 0.00
#&gt; SP101666     3  0.3983     0.5341 0.00 0.00 0.66 0.34 0.00
#&gt; SP101670     2  0.4307     0.0559 0.00 0.50 0.00 0.00 0.50
#&gt; SP101674     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101678     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101682     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101686     2  0.3561     0.6449 0.00 0.74 0.00 0.00 0.26
#&gt; SP101690     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101694     2  0.3895     0.5586 0.00 0.68 0.00 0.00 0.32
#&gt; SP101700     3  0.3291     0.7287 0.00 0.00 0.84 0.04 0.12
#&gt; SP101708     3  0.0609     0.7690 0.00 0.00 0.98 0.02 0.00
#&gt; SP101716     2  0.4060     0.4773 0.00 0.64 0.36 0.00 0.00
#&gt; SP101724     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101732     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101740     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101795     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101845     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101881     3  0.0609     0.7690 0.00 0.00 0.98 0.02 0.00
#&gt; SP101891     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP101921     5  0.0000     0.8625 0.00 0.00 0.00 0.00 1.00
#&gt; SP101931     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102015     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102035     3  0.2280     0.7348 0.00 0.00 0.88 0.12 0.00
#&gt; SP102045     2  0.2280     0.8116 0.00 0.88 0.00 0.00 0.12
#&gt; SP102055     3  0.4126     0.4189 0.00 0.00 0.62 0.38 0.00
#&gt; SP102064     3  0.3561     0.6052 0.00 0.00 0.74 0.26 0.00
#&gt; SP102074     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102084     3  0.0609     0.7690 0.00 0.00 0.98 0.02 0.00
#&gt; SP102090     2  0.3274     0.7053 0.00 0.78 0.22 0.00 0.00
#&gt; SP102096     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102103     2  0.2280     0.8095 0.00 0.88 0.00 0.00 0.12
#&gt; SP102113     3  0.0609     0.7690 0.00 0.00 0.98 0.02 0.00
#&gt; SP102123     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102133     2  0.2516     0.7931 0.00 0.86 0.00 0.00 0.14
#&gt; SP102143     3  0.4126     0.4778 0.00 0.00 0.62 0.38 0.00
#&gt; SP102161     5  0.0000     0.8625 0.00 0.00 0.00 0.00 1.00
#&gt; SP102168     3  0.3109     0.6372 0.00 0.20 0.80 0.00 0.00
#&gt; SP102174     2  0.0000     0.8742 0.00 1.00 0.00 0.00 0.00
#&gt; SP102187     3  0.0609     0.7715 0.02 0.00 0.98 0.00 0.00
#&gt; SP102485     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102489     1  0.3895     0.5339 0.68 0.00 0.00 0.00 0.32
#&gt; SP102499     1  0.4458     0.6878 0.76 0.00 0.12 0.12 0.00
#&gt; SP102507     4  0.0609     0.7974 0.02 0.00 0.00 0.98 0.00
#&gt; SP102511     1  0.5032     0.6630 0.74 0.00 0.12 0.12 0.02
#&gt; SP102517     5  0.2929     0.7347 0.18 0.00 0.00 0.00 0.82
#&gt; SP102523     5  0.3274     0.6643 0.22 0.00 0.00 0.00 0.78
#&gt; SP102529     5  0.1410     0.8548 0.06 0.00 0.00 0.00 0.94
#&gt; SP102537     3  0.5599     0.5549 0.26 0.00 0.62 0.12 0.00
#&gt; SP102541     5  0.0000     0.8625 0.00 0.00 0.00 0.00 1.00
#&gt; SP102547     1  0.4840     0.4765 0.64 0.00 0.00 0.04 0.32
#&gt; SP102557     1  0.4794     0.6938 0.76 0.00 0.10 0.12 0.02
#&gt; SP102561     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102567     5  0.1410     0.8548 0.06 0.00 0.00 0.00 0.94
#&gt; SP102573     5  0.3274     0.6643 0.22 0.00 0.00 0.00 0.78
#&gt; SP102581     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102591     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102597     1  0.4458     0.6878 0.76 0.00 0.12 0.12 0.00
#&gt; SP102605     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102611     5  0.3684     0.5782 0.00 0.28 0.00 0.00 0.72
#&gt; SP102617     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102620     4  0.0000     0.7986 0.00 0.00 0.00 1.00 0.00
#&gt; SP102622     5  0.1410     0.8548 0.06 0.00 0.00 0.00 0.94
#&gt; SP102626     4  0.0000     0.7986 0.00 0.00 0.00 1.00 0.00
#&gt; SP102630     4  0.2516     0.7219 0.00 0.00 0.14 0.86 0.00
#&gt; SP102633     4  0.1043     0.8048 0.04 0.00 0.00 0.96 0.00
#&gt; SP102647     4  0.0000     0.7986 0.00 0.00 0.00 1.00 0.00
#&gt; SP102652     4  0.3513     0.7064 0.00 0.00 0.02 0.80 0.18
#&gt; SP102690     1  0.0609     0.9351 0.98 0.00 0.00 0.02 0.00
#&gt; SP102716     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102718     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102733     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102741     1  0.1043     0.9189 0.96 0.00 0.00 0.04 0.00
#&gt; SP102747     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102755     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102759     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102783     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102804     4  0.3109     0.7488 0.20 0.00 0.00 0.80 0.00
#&gt; SP102816     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102827     5  0.3684     0.6042 0.28 0.00 0.00 0.00 0.72
#&gt; SP102839     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102873     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102881     1  0.1043     0.9183 0.96 0.00 0.00 0.04 0.00
#&gt; SP102897     1  0.1043     0.9182 0.96 0.00 0.00 0.04 0.00
#&gt; SP102913     5  0.0000     0.8625 0.00 0.00 0.00 0.00 1.00
#&gt; SP102921     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102929     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102945     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102957     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt; SP102965     1  0.0000     0.9524 1.00 0.00 0.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2612 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-4-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-5'>
<p><a id='tab-ATC-skmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5   p6
#&gt; SP1003       2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP10084      2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP1009       4  0.0000     0.6108 0.00 0.00 0.00 1.00 0.00 0.00
#&gt; SP10150      2  0.4282     0.3508 0.00 0.56 0.00 0.00 0.42 0.02
#&gt; SP101515     3  0.5012     0.5644 0.00 0.00 0.60 0.30 0.00 0.10
#&gt; SP101519     5  0.5260     0.1120 0.00 0.44 0.02 0.02 0.50 0.02
#&gt; SP101521     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101523     3  0.2094     0.7372 0.02 0.00 0.90 0.00 0.00 0.08
#&gt; SP101526     2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP101528     4  0.5747     0.5461 0.00 0.00 0.20 0.50 0.00 0.30
#&gt; SP101532     4  0.5992    -0.1649 0.02 0.00 0.14 0.48 0.00 0.36
#&gt; SP101536     2  0.4310     0.2816 0.00 0.54 0.00 0.00 0.44 0.02
#&gt; SP101540     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101544     2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP101548     4  0.1267     0.6111 0.00 0.00 0.00 0.94 0.00 0.06
#&gt; SP101552     2  0.0547     0.8683 0.00 0.98 0.00 0.00 0.02 0.00
#&gt; SP101558     3  0.2094     0.7372 0.02 0.00 0.90 0.00 0.00 0.08
#&gt; SP101564     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101572     2  0.3163     0.7910 0.00 0.82 0.00 0.00 0.14 0.04
#&gt; SP101576     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101580     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101584     2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP101588     2  0.5262     0.4770 0.00 0.62 0.04 0.02 0.30 0.02
#&gt; SP101592     2  0.0547     0.8684 0.00 0.98 0.00 0.00 0.00 0.02
#&gt; SP101596     2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP101600     3  0.2094     0.7372 0.02 0.00 0.90 0.00 0.00 0.08
#&gt; SP101604     4  0.5371     0.5128 0.00 0.00 0.12 0.52 0.00 0.36
#&gt; SP101610     6  0.5882    -0.1265 0.00 0.00 0.38 0.20 0.00 0.42
#&gt; SP101616     3  0.3797     0.3028 0.00 0.42 0.58 0.00 0.00 0.00
#&gt; SP101622     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101628     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101634     3  0.5888     0.3398 0.00 0.00 0.40 0.40 0.00 0.20
#&gt; SP101642     2  0.1814     0.8088 0.00 0.90 0.10 0.00 0.00 0.00
#&gt; SP101648     3  0.2350     0.7531 0.00 0.00 0.88 0.02 0.00 0.10
#&gt; SP101654     2  0.0937     0.8558 0.00 0.96 0.04 0.00 0.00 0.00
#&gt; SP101658     6  0.5529     0.6698 0.20 0.00 0.12 0.00 0.04 0.64
#&gt; SP101662     6  0.5324    -0.1935 0.00 0.00 0.12 0.34 0.00 0.54
#&gt; SP101666     6  0.5876     0.2943 0.00 0.00 0.26 0.26 0.00 0.48
#&gt; SP101670     2  0.4310     0.1809 0.00 0.54 0.00 0.00 0.44 0.02
#&gt; SP101674     2  0.0547     0.8684 0.00 0.98 0.00 0.00 0.00 0.02
#&gt; SP101678     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101682     2  0.0937     0.8667 0.00 0.96 0.00 0.00 0.00 0.04
#&gt; SP101686     2  0.2793     0.7303 0.00 0.80 0.00 0.00 0.20 0.00
#&gt; SP101690     2  0.0547     0.8652 0.00 0.98 0.02 0.00 0.00 0.00
#&gt; SP101694     2  0.3460     0.7148 0.00 0.76 0.00 0.00 0.22 0.02
#&gt; SP101700     6  0.6200     0.1394 0.00 0.00 0.34 0.18 0.02 0.46
#&gt; SP101708     3  0.1814     0.7433 0.00 0.00 0.90 0.00 0.00 0.10
#&gt; SP101716     3  0.3828     0.2426 0.00 0.44 0.56 0.00 0.00 0.00
#&gt; SP101724     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101732     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101740     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101795     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101845     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101881     3  0.2474     0.7546 0.00 0.00 0.88 0.04 0.00 0.08
#&gt; SP101891     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP101921     5  0.2725     0.8082 0.00 0.00 0.02 0.04 0.88 0.06
#&gt; SP101931     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102015     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102035     3  0.4873     0.6040 0.00 0.00 0.60 0.32 0.00 0.08
#&gt; SP102045     2  0.2790     0.7961 0.00 0.84 0.00 0.00 0.14 0.02
#&gt; SP102055     3  0.5071     0.4294 0.00 0.00 0.52 0.40 0.00 0.08
#&gt; SP102064     6  0.3660     0.4450 0.00 0.00 0.16 0.06 0.00 0.78
#&gt; SP102074     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102084     3  0.1814     0.7433 0.00 0.00 0.90 0.00 0.00 0.10
#&gt; SP102090     3  0.3851     0.1749 0.00 0.46 0.54 0.00 0.00 0.00
#&gt; SP102096     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102103     2  0.1814     0.8266 0.00 0.90 0.00 0.00 0.10 0.00
#&gt; SP102113     3  0.1814     0.7433 0.00 0.00 0.90 0.00 0.00 0.10
#&gt; SP102123     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102133     2  0.1814     0.8266 0.00 0.90 0.00 0.00 0.10 0.00
#&gt; SP102143     4  0.5371    -0.0521 0.00 0.00 0.12 0.52 0.00 0.36
#&gt; SP102161     5  0.0000     0.8358 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; SP102168     3  0.3163     0.7235 0.00 0.14 0.82 0.00 0.00 0.04
#&gt; SP102174     2  0.0000     0.8689 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; SP102187     3  0.2094     0.7372 0.02 0.00 0.90 0.00 0.00 0.08
#&gt; SP102485     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102489     1  0.5945     0.2195 0.56 0.00 0.02 0.04 0.32 0.06
#&gt; SP102499     6  0.4798     0.6855 0.30 0.00 0.08 0.00 0.00 0.62
#&gt; SP102507     4  0.0937     0.6140 0.00 0.00 0.00 0.96 0.00 0.04
#&gt; SP102511     6  0.4798     0.6855 0.30 0.00 0.08 0.00 0.00 0.62
#&gt; SP102517     5  0.2793     0.6600 0.20 0.00 0.00 0.00 0.80 0.00
#&gt; SP102523     5  0.5214     0.6248 0.18 0.00 0.02 0.04 0.70 0.06
#&gt; SP102529     5  0.1267     0.8265 0.06 0.00 0.00 0.00 0.94 0.00
#&gt; SP102537     6  0.5029     0.6958 0.26 0.00 0.12 0.00 0.00 0.62
#&gt; SP102541     5  0.2725     0.8082 0.00 0.00 0.02 0.04 0.88 0.06
#&gt; SP102547     1  0.6242     0.1948 0.56 0.00 0.02 0.04 0.28 0.10
#&gt; SP102557     6  0.4798     0.6855 0.30 0.00 0.08 0.00 0.00 0.62
#&gt; SP102561     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102567     5  0.1267     0.8265 0.06 0.00 0.00 0.00 0.94 0.00
#&gt; SP102573     5  0.5214     0.6248 0.18 0.00 0.02 0.04 0.70 0.06
#&gt; SP102581     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102591     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102597     6  0.4798     0.6855 0.30 0.00 0.08 0.00 0.00 0.62
#&gt; SP102605     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102611     5  0.3198     0.6227 0.00 0.26 0.00 0.00 0.74 0.00
#&gt; SP102617     1  0.0937     0.8876 0.96 0.00 0.00 0.00 0.00 0.04
#&gt; SP102620     4  0.2048     0.6280 0.00 0.00 0.00 0.88 0.00 0.12
#&gt; SP102622     5  0.1267     0.8265 0.06 0.00 0.00 0.00 0.94 0.00
#&gt; SP102626     4  0.0937     0.6140 0.00 0.00 0.00 0.96 0.00 0.04
#&gt; SP102630     4  0.2631     0.4262 0.00 0.00 0.18 0.82 0.00 0.00
#&gt; SP102633     4  0.2474     0.6322 0.04 0.00 0.00 0.88 0.00 0.08
#&gt; SP102647     4  0.0937     0.6140 0.00 0.00 0.00 0.96 0.00 0.04
#&gt; SP102652     4  0.3942     0.5043 0.00 0.00 0.04 0.80 0.10 0.06
#&gt; SP102690     1  0.3711     0.5573 0.72 0.00 0.00 0.02 0.00 0.26
#&gt; SP102716     1  0.0937     0.8875 0.96 0.00 0.00 0.00 0.00 0.04
#&gt; SP102718     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102733     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102741     1  0.2631     0.6932 0.82 0.00 0.00 0.00 0.00 0.18
#&gt; SP102747     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102755     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102759     1  0.0000     0.9277 1.00 0.00 0.00 0.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2626 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-5-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-6'>
<p><a id='tab-ATC-skmeans-get-classes-6-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 7), get_membership(res, k = 7))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5   p6   p7
#&gt; SP1003       2  0.3199     0.7894 0.00 0.80 0.00 0.00 0.06 0.00 0.14
#&gt; SP10084      2  0.2569     0.8027 0.00 0.84 0.00 0.00 0.02 0.00 0.14
#&gt; SP1009       6  0.6264     0.2179 0.00 0.00 0.04 0.30 0.00 0.36 0.30
#&gt; SP10150      2  0.4945     0.3786 0.00 0.52 0.00 0.00 0.36 0.00 0.12
#&gt; SP101515     3  0.5086     0.5469 0.00 0.00 0.64 0.06 0.00 0.08 0.22
#&gt; SP101519     7  0.5332     0.2424 0.00 0.38 0.00 0.00 0.18 0.00 0.44
#&gt; SP101521     2  0.0504     0.8168 0.00 0.98 0.00 0.00 0.00 0.00 0.02
#&gt; SP101523     3  0.2569     0.7303 0.02 0.00 0.84 0.00 0.00 0.14 0.00
#&gt; SP101526     2  0.2376     0.8085 0.00 0.86 0.00 0.00 0.02 0.00 0.12
#&gt; SP101528     4  0.4535     0.4653 0.00 0.00 0.24 0.64 0.00 0.00 0.12
#&gt; SP101532     6  0.5840     0.3050 0.02 0.00 0.12 0.04 0.00 0.56 0.26
#&gt; SP101536     5  0.4487     0.0556 0.00 0.42 0.00 0.00 0.52 0.00 0.06
#&gt; SP101540     2  0.0863     0.8168 0.00 0.96 0.00 0.00 0.00 0.00 0.04
#&gt; SP101544     2  0.2569     0.8027 0.00 0.84 0.00 0.00 0.02 0.00 0.14
#&gt; SP101548     6  0.6413     0.2310 0.00 0.00 0.06 0.30 0.00 0.38 0.26
#&gt; SP101552     2  0.0504     0.8136 0.00 0.98 0.00 0.00 0.00 0.00 0.02
#&gt; SP101558     3  0.2569     0.7303 0.02 0.00 0.84 0.00 0.00 0.14 0.00
#&gt; SP101564     2  0.0504     0.8157 0.00 0.98 0.00 0.00 0.02 0.00 0.00
#&gt; SP101572     2  0.4478     0.6811 0.00 0.66 0.00 0.00 0.20 0.00 0.14
#&gt; SP101576     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101580     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101584     2  0.2569     0.8027 0.00 0.84 0.00 0.00 0.02 0.00 0.14
#&gt; SP101588     2  0.4264     0.4397 0.00 0.62 0.00 0.00 0.06 0.00 0.32
#&gt; SP101592     2  0.1664     0.8151 0.00 0.92 0.00 0.00 0.02 0.00 0.06
#&gt; SP101596     2  0.2376     0.8083 0.00 0.86 0.00 0.00 0.02 0.00 0.12
#&gt; SP101600     3  0.2569     0.7303 0.02 0.00 0.84 0.00 0.00 0.14 0.00
#&gt; SP101604     4  0.5681     0.2548 0.00 0.00 0.08 0.58 0.00 0.12 0.22
#&gt; SP101610     3  0.6952     0.1398 0.00 0.00 0.34 0.26 0.00 0.24 0.16
#&gt; SP101616     2  0.3562    -0.0912 0.00 0.50 0.50 0.00 0.00 0.00 0.00
#&gt; SP101622     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101628     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101634     6  0.5604     0.1524 0.00 0.00 0.26 0.02 0.00 0.52 0.20
#&gt; SP101642     2  0.1671     0.7714 0.00 0.90 0.10 0.00 0.00 0.00 0.00
#&gt; SP101648     3  0.2569     0.7343 0.00 0.00 0.84 0.02 0.00 0.14 0.00
#&gt; SP101654     2  0.0863     0.8085 0.00 0.96 0.04 0.00 0.00 0.00 0.00
#&gt; SP101658     6  0.5363     0.3682 0.22 0.00 0.00 0.14 0.00 0.60 0.04
#&gt; SP101662     4  0.5317     0.1405 0.00 0.00 0.02 0.56 0.00 0.28 0.14
#&gt; SP101666     6  0.3388     0.1428 0.00 0.00 0.20 0.04 0.00 0.76 0.00
#&gt; SP101670     2  0.4970     0.3235 0.00 0.58 0.00 0.00 0.18 0.00 0.24
#&gt; SP101674     2  0.1664     0.8151 0.00 0.92 0.00 0.00 0.02 0.00 0.06
#&gt; SP101678     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101682     2  0.2376     0.8067 0.00 0.86 0.00 0.00 0.02 0.00 0.12
#&gt; SP101686     2  0.2376     0.7521 0.00 0.86 0.00 0.00 0.12 0.00 0.02
#&gt; SP101690     2  0.0504     0.8146 0.00 0.98 0.02 0.00 0.00 0.00 0.00
#&gt; SP101694     2  0.4015     0.6395 0.00 0.68 0.00 0.00 0.26 0.00 0.06
#&gt; SP101700     6  0.5259    -0.0152 0.00 0.00 0.22 0.00 0.00 0.52 0.26
#&gt; SP101708     3  0.2569     0.7343 0.00 0.00 0.84 0.02 0.00 0.14 0.00
#&gt; SP101716     3  0.3546     0.1866 0.00 0.46 0.54 0.00 0.00 0.00 0.00
#&gt; SP101724     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101732     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101740     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101795     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101845     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101881     3  0.2259     0.7354 0.00 0.00 0.84 0.00 0.00 0.16 0.00
#&gt; SP101891     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP101921     7  0.3459     0.6177 0.00 0.00 0.00 0.00 0.40 0.00 0.60
#&gt; SP101931     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102015     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102035     3  0.5452     0.2084 0.00 0.00 0.46 0.00 0.00 0.30 0.24
#&gt; SP102045     2  0.4015     0.6319 0.00 0.68 0.00 0.00 0.26 0.00 0.06
#&gt; SP102055     3  0.6215     0.2122 0.00 0.00 0.46 0.06 0.00 0.26 0.22
#&gt; SP102064     6  0.3867     0.1982 0.00 0.00 0.02 0.38 0.00 0.60 0.00
#&gt; SP102074     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102084     3  0.2569     0.7343 0.00 0.00 0.84 0.02 0.00 0.14 0.00
#&gt; SP102090     2  0.3562    -0.0654 0.00 0.50 0.50 0.00 0.00 0.00 0.00
#&gt; SP102096     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102103     2  0.1433     0.7939 0.00 0.92 0.00 0.00 0.08 0.00 0.00
#&gt; SP102113     3  0.2569     0.7349 0.00 0.00 0.84 0.00 0.00 0.14 0.02
#&gt; SP102123     2  0.0000     0.8157 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102133     2  0.1886     0.7721 0.00 0.88 0.00 0.00 0.12 0.00 0.00
#&gt; SP102143     6  0.5416     0.1632 0.00 0.00 0.20 0.00 0.00 0.42 0.38
#&gt; SP102161     5  0.1166     0.7217 0.00 0.00 0.00 0.00 0.94 0.00 0.06
#&gt; SP102168     3  0.3449     0.7064 0.00 0.14 0.78 0.00 0.00 0.08 0.00
#&gt; SP102174     2  0.0504     0.8188 0.00 0.98 0.00 0.00 0.00 0.00 0.02
#&gt; SP102187     3  0.2569     0.7303 0.02 0.00 0.84 0.00 0.00 0.14 0.00
#&gt; SP102485     1  0.0000     0.9525 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102489     7  0.4873     0.5191 0.18 0.00 0.00 0.00 0.22 0.00 0.60
#&gt; SP102499     6  0.5055     0.3741 0.26 0.00 0.00 0.18 0.00 0.56 0.00
#&gt; SP102507     6  0.6627     0.2253 0.02 0.00 0.04 0.32 0.00 0.36 0.26
#&gt; SP102511     6  0.5417     0.3750 0.24 0.00 0.00 0.18 0.02 0.56 0.00
#&gt; SP102517     5  0.1433     0.7379 0.08 0.00 0.00 0.00 0.92 0.00 0.00
#&gt; SP102523     7  0.4127     0.6323 0.04 0.00 0.00 0.00 0.36 0.00 0.60
#&gt; SP102529     5  0.1166     0.7678 0.06 0.00 0.00 0.00 0.94 0.00 0.00
#&gt; SP102537     6  0.4789     0.3727 0.26 0.00 0.00 0.14 0.00 0.60 0.00
#&gt; SP102541     7  0.3496     0.5994 0.00 0.00 0.00 0.00 0.42 0.00 0.58
#&gt; SP102547     1  0.6020     0.1998 0.54 0.00 0.00 0.10 0.22 0.00 0.14
#&gt; SP102557     6  0.5417     0.3750 0.24 0.00 0.00 0.18 0.02 0.56 0.00
#&gt; SP102561     1  0.0000     0.9525 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102567     5  0.1166     0.7678 0.06 0.00 0.00 0.00 0.94 0.00 0.00
#&gt; SP102573     7  0.4127     0.6323 0.04 0.00 0.00 0.00 0.36 0.00 0.60
#&gt; SP102581     1  0.0000     0.9525 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102591     1  0.0000     0.9525 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102597     6  0.5055     0.3741 0.26 0.00 0.00 0.18 0.00 0.56 0.00
#&gt; SP102605     1  0.0000     0.9525 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102611     5  0.3139     0.5043 0.00 0.30 0.00 0.00 0.70 0.00 0.00
#&gt; SP102617     1  0.2376     0.8109 0.86 0.00 0.00 0.12 0.00 0.02 0.00
#&gt; SP102620     4  0.5579    -0.1574 0.00 0.00 0.00 0.38 0.00 0.36 0.26
#&gt; SP102622     5  0.1166     0.7678 0.06 0.00 0.00 0.00 0.94 0.00 0.00
#&gt; SP102626     6  0.6237     0.2161 0.00 0.00 0.04 0.34 0.00 0.36 0.26
#&gt; SP102630     6  0.6752     0.2480 0.00 0.00 0.20 0.14 0.00 0.40 0.26
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2637 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-6-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-6-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-6-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-7'>
<p><a id='tab-ATC-skmeans-get-classes-7-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 8), get_membership(res, k = 8))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5   p6   p7   p8
#&gt; SP1003       2   0.286     0.6377 0.00 0.84 0.00 0.06 0.08 0.02 0.00 0.00
#&gt; SP10084      2   0.141     0.6563 0.00 0.94 0.00 0.02 0.02 0.02 0.00 0.00
#&gt; SP1009       8   0.128     0.7354 0.00 0.00 0.00 0.02 0.00 0.00 0.04 0.94
#&gt; SP10150      2   0.467     0.4517 0.00 0.60 0.00 0.00 0.28 0.02 0.10 0.00
#&gt; SP101515     3   0.512     0.3216 0.00 0.00 0.56 0.04 0.00 0.06 0.02 0.32
#&gt; SP101519     7   0.504     0.4916 0.00 0.12 0.02 0.06 0.06 0.04 0.70 0.00
#&gt; SP101521     2   0.434     0.6867 0.00 0.70 0.00 0.16 0.00 0.04 0.10 0.00
#&gt; SP101523     3   0.272     0.6416 0.02 0.00 0.80 0.00 0.00 0.18 0.00 0.00
#&gt; SP101526     2   0.174     0.6765 0.00 0.92 0.00 0.04 0.02 0.02 0.00 0.00
#&gt; SP101528     4   0.573     0.1823 0.00 0.00 0.22 0.42 0.00 0.04 0.00 0.32
#&gt; SP101532     8   0.387     0.6160 0.00 0.00 0.06 0.02 0.00 0.20 0.00 0.72
#&gt; SP101536     5   0.357     0.2787 0.00 0.36 0.00 0.00 0.62 0.02 0.00 0.00
#&gt; SP101540     2   0.505     0.6890 0.00 0.66 0.02 0.14 0.00 0.06 0.12 0.00
#&gt; SP101544     2   0.262     0.6486 0.00 0.86 0.00 0.06 0.02 0.06 0.00 0.00
#&gt; SP101548     8   0.259     0.7157 0.00 0.00 0.04 0.02 0.00 0.08 0.00 0.86
#&gt; SP101552     2   0.577     0.6695 0.00 0.60 0.02 0.14 0.02 0.06 0.16 0.00
#&gt; SP101558     3   0.272     0.6416 0.02 0.00 0.80 0.00 0.00 0.18 0.00 0.00
#&gt; SP101564     2   0.483     0.6948 0.00 0.72 0.02 0.10 0.04 0.04 0.08 0.00
#&gt; SP101572     2   0.317     0.5809 0.00 0.78 0.00 0.02 0.18 0.02 0.00 0.00
#&gt; SP101576     2   0.452     0.6937 0.00 0.74 0.02 0.10 0.02 0.04 0.08 0.00
#&gt; SP101580     2   0.506     0.6876 0.00 0.68 0.02 0.14 0.02 0.04 0.10 0.00
#&gt; SP101584     2   0.156     0.6780 0.00 0.92 0.00 0.06 0.02 0.00 0.00 0.00
#&gt; SP101588     7   0.629    -0.0418 0.00 0.24 0.02 0.16 0.02 0.06 0.50 0.00
#&gt; SP101592     2   0.305     0.6807 0.00 0.84 0.00 0.04 0.08 0.02 0.02 0.00
#&gt; SP101596     2   0.128     0.6828 0.00 0.94 0.00 0.00 0.02 0.00 0.04 0.00
#&gt; SP101600     3   0.272     0.6416 0.02 0.00 0.80 0.00 0.00 0.18 0.00 0.00
#&gt; SP101604     8   0.550     0.2570 0.00 0.00 0.06 0.22 0.00 0.18 0.00 0.54
#&gt; SP101610     8   0.654     0.1285 0.00 0.00 0.26 0.20 0.00 0.20 0.00 0.34
#&gt; SP101616     3   0.659    -0.2462 0.00 0.30 0.38 0.18 0.00 0.02 0.12 0.00
#&gt; SP101622     2   0.549     0.6749 0.00 0.64 0.02 0.14 0.02 0.06 0.12 0.00
#&gt; SP101628     2   0.549     0.6749 0.00 0.64 0.02 0.14 0.02 0.06 0.12 0.00
#&gt; SP101634     8   0.490     0.3674 0.00 0.00 0.28 0.00 0.00 0.20 0.00 0.52
#&gt; SP101642     2   0.683     0.6195 0.00 0.50 0.10 0.18 0.02 0.08 0.12 0.00
#&gt; SP101648     3   0.286     0.6627 0.00 0.00 0.82 0.02 0.00 0.14 0.00 0.02
#&gt; SP101654     2   0.591     0.6668 0.00 0.60 0.04 0.16 0.02 0.06 0.12 0.00
#&gt; SP101658     6   0.428     0.6578 0.18 0.00 0.00 0.00 0.00 0.70 0.08 0.04
#&gt; SP101662     6   0.623    -0.0516 0.00 0.00 0.10 0.22 0.00 0.36 0.00 0.32
#&gt; SP101666     6   0.461     0.2746 0.00 0.00 0.16 0.00 0.00 0.58 0.00 0.26
#&gt; SP101670     7   0.677    -0.0816 0.00 0.24 0.02 0.16 0.06 0.06 0.46 0.00
#&gt; SP101674     2   0.422     0.6934 0.00 0.76 0.00 0.08 0.06 0.04 0.06 0.00
#&gt; SP101678     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101682     2   0.128     0.6683 0.00 0.94 0.00 0.00 0.04 0.02 0.00 0.00
#&gt; SP101686     2   0.671     0.5733 0.00 0.48 0.02 0.16 0.06 0.06 0.22 0.00
#&gt; SP101690     2   0.531     0.6812 0.00 0.66 0.02 0.14 0.02 0.06 0.10 0.00
#&gt; SP101694     2   0.294     0.5139 0.00 0.70 0.00 0.00 0.30 0.00 0.00 0.00
#&gt; SP101700     6   0.527     0.2858 0.00 0.00 0.18 0.00 0.00 0.56 0.22 0.04
#&gt; SP101708     3   0.286     0.6627 0.00 0.00 0.82 0.02 0.00 0.14 0.00 0.02
#&gt; SP101716     3   0.660    -0.0514 0.00 0.24 0.44 0.18 0.00 0.04 0.10 0.00
#&gt; SP101724     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101732     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101740     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101795     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101845     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101881     3   0.273     0.6642 0.00 0.00 0.82 0.00 0.00 0.14 0.00 0.04
#&gt; SP101891     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP101921     7   0.195     0.7073 0.00 0.00 0.00 0.00 0.14 0.00 0.86 0.00
#&gt; SP101931     2   0.332     0.7033 0.00 0.82 0.00 0.08 0.02 0.02 0.06 0.00
#&gt; SP102015     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP102035     8   0.451     0.2145 0.00 0.00 0.38 0.00 0.00 0.10 0.00 0.52
#&gt; SP102045     2   0.351     0.4557 0.00 0.64 0.00 0.00 0.34 0.02 0.00 0.00
#&gt; SP102055     8   0.507     0.0357 0.00 0.00 0.44 0.04 0.00 0.08 0.00 0.44
#&gt; SP102064     6   0.476     0.4204 0.00 0.00 0.08 0.22 0.00 0.64 0.00 0.06
#&gt; SP102074     2   0.562     0.6697 0.00 0.62 0.02 0.16 0.02 0.06 0.12 0.00
#&gt; SP102084     3   0.272     0.6503 0.00 0.00 0.80 0.02 0.00 0.18 0.00 0.00
#&gt; SP102090     3   0.735    -0.1413 0.00 0.28 0.36 0.16 0.02 0.10 0.08 0.00
#&gt; SP102096     2   0.490     0.6880 0.00 0.70 0.02 0.12 0.02 0.04 0.10 0.00
#&gt; SP102103     2   0.630     0.6501 0.00 0.56 0.02 0.16 0.06 0.06 0.14 0.00
#&gt; SP102113     3   0.286     0.6636 0.00 0.00 0.82 0.00 0.00 0.14 0.02 0.02
#&gt; SP102123     2   0.577     0.6688 0.00 0.60 0.02 0.16 0.02 0.06 0.14 0.00
#&gt; SP102133     2   0.681     0.6178 0.00 0.50 0.02 0.16 0.14 0.06 0.12 0.00
#&gt; SP102143     8   0.536     0.4850 0.00 0.00 0.14 0.00 0.00 0.20 0.08 0.58
#&gt; SP102161     5   0.156     0.7327 0.00 0.00 0.00 0.00 0.90 0.00 0.10 0.00
#&gt; SP102168     3   0.395     0.6366 0.00 0.08 0.76 0.02 0.00 0.12 0.02 0.00
#&gt; SP102174     2   0.452     0.6998 0.00 0.74 0.02 0.10 0.02 0.04 0.08 0.00
#&gt; SP102187     3   0.272     0.6416 0.02 0.00 0.80 0.00 0.00 0.18 0.00 0.00
#&gt; SP102485     1   0.000     0.9360 1.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102489     7   0.241     0.6812 0.08 0.00 0.00 0.00 0.06 0.00 0.86 0.00
#&gt; SP102499     6   0.470     0.6828 0.26 0.00 0.00 0.08 0.00 0.62 0.00 0.04
#&gt; SP102507     8   0.134     0.7210 0.00 0.00 0.00 0.08 0.00 0.00 0.00 0.92
#&gt; SP102511     6   0.524     0.6899 0.22 0.00 0.00 0.08 0.04 0.62 0.00 0.04
#&gt; SP102517     5   0.156     0.7710 0.06 0.00 0.00 0.00 0.92 0.00 0.02 0.00
#&gt; SP102523     7   0.195     0.7159 0.00 0.00 0.00 0.00 0.14 0.00 0.86 0.00
#&gt; SP102529     5   0.156     0.7710 0.06 0.00 0.00 0.00 0.92 0.00 0.02 0.00
#&gt; SP102537     6   0.376     0.6903 0.26 0.00 0.00 0.06 0.00 0.68 0.00 0.00
#&gt; SP102541     7   0.211     0.7054 0.00 0.00 0.00 0.00 0.16 0.00 0.84 0.00
#&gt; SP102547     1   0.575    -0.1648 0.42 0.00 0.00 0.06 0.06 0.04 0.42 0.00
#&gt; SP102557     6   0.524     0.6899 0.22 0.00 0.00 0.08 0.04 0.62 0.00 0.04
#&gt; SP102561     1   0.000     0.9360 1.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102567     5   0.156     0.7710 0.06 0.00 0.00 0.00 0.92 0.00 0.02 0.00
#&gt; SP102573     7   0.195     0.7159 0.00 0.00 0.00 0.00 0.14 0.00 0.86 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2647 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-7-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-7-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-7-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-8'>
<p><a id='tab-ATC-skmeans-get-classes-8-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 9), get_membership(res, k = 9))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5   p6   p7   p8   p9
#&gt; SP1003       2  0.1473     0.6274 0.00 0.92 0.00 0.00 0.02 0.00 0.00 0.00 0.06
#&gt; SP10084      2  0.0446     0.6327 0.00 0.98 0.00 0.00 0.00 0.00 0.00 0.00 0.02
#&gt; SP1009       8  0.1786     0.7280 0.00 0.00 0.00 0.06 0.00 0.00 0.04 0.90 0.00
#&gt; SP10150      2  0.3758     0.5095 0.00 0.68 0.00 0.00 0.22 0.00 0.10 0.00 0.00
#&gt; SP101515     3  0.5397     0.3803 0.00 0.00 0.58 0.06 0.00 0.08 0.02 0.24 0.02
#&gt; SP101519     7  0.4758     0.4469 0.00 0.06 0.00 0.00 0.06 0.00 0.52 0.00 0.36
#&gt; SP101521     9  0.3063     0.4113 0.00 0.40 0.00 0.00 0.00 0.00 0.00 0.00 0.60
#&gt; SP101523     3  0.2431     0.6011 0.02 0.00 0.82 0.00 0.00 0.16 0.00 0.00 0.00
#&gt; SP101526     2  0.1916     0.5898 0.00 0.88 0.00 0.00 0.02 0.00 0.00 0.00 0.10
#&gt; SP101528     8  0.6764    -0.0534 0.00 0.00 0.14 0.34 0.00 0.10 0.06 0.34 0.02
#&gt; SP101532     8  0.3637     0.5788 0.00 0.00 0.02 0.04 0.00 0.24 0.00 0.70 0.00
#&gt; SP101536     5  0.4110     0.3805 0.00 0.18 0.00 0.00 0.64 0.00 0.00 0.00 0.18
#&gt; SP101540     2  0.3140    -0.1669 0.00 0.54 0.00 0.00 0.00 0.00 0.00 0.00 0.46
#&gt; SP101544     2  0.1670     0.5989 0.00 0.88 0.00 0.00 0.00 0.00 0.00 0.00 0.12
#&gt; SP101548     8  0.4619     0.6349 0.00 0.00 0.02 0.10 0.00 0.14 0.06 0.68 0.00
#&gt; SP101552     9  0.3424     0.4897 0.00 0.38 0.00 0.00 0.02 0.00 0.00 0.00 0.60
#&gt; SP101558     3  0.2275     0.6086 0.02 0.00 0.84 0.00 0.00 0.14 0.00 0.00 0.00
#&gt; SP101564     2  0.3997    -0.2230 0.00 0.48 0.00 0.00 0.06 0.00 0.00 0.00 0.46
#&gt; SP101572     2  0.1843     0.6131 0.00 0.86 0.00 0.00 0.14 0.00 0.00 0.00 0.00
#&gt; SP101576     2  0.3990    -0.1727 0.00 0.50 0.00 0.00 0.06 0.00 0.00 0.00 0.44
#&gt; SP101580     9  0.3778     0.3858 0.00 0.44 0.00 0.00 0.04 0.00 0.00 0.00 0.52
#&gt; SP101584     2  0.2001     0.5823 0.00 0.84 0.00 0.00 0.00 0.00 0.00 0.00 0.16
#&gt; SP101588     9  0.2703     0.5220 0.00 0.02 0.00 0.00 0.00 0.00 0.20 0.00 0.78
#&gt; SP101592     2  0.3122     0.4944 0.00 0.76 0.00 0.00 0.06 0.00 0.00 0.00 0.18
#&gt; SP101596     2  0.2224     0.5888 0.00 0.86 0.00 0.00 0.04 0.00 0.00 0.00 0.10
#&gt; SP101600     3  0.2275     0.6086 0.02 0.00 0.84 0.00 0.00 0.14 0.00 0.00 0.00
#&gt; SP101604     8  0.6727     0.2359 0.00 0.00 0.14 0.22 0.00 0.22 0.06 0.36 0.00
#&gt; SP101610     8  0.6833     0.2086 0.00 0.00 0.20 0.20 0.00 0.20 0.06 0.34 0.00
#&gt; SP101616     9  0.3572     0.5514 0.00 0.08 0.22 0.00 0.00 0.00 0.00 0.00 0.70
#&gt; SP101622     9  0.3478     0.5909 0.00 0.30 0.00 0.00 0.04 0.00 0.00 0.00 0.66
#&gt; SP101628     9  0.3604     0.5521 0.00 0.34 0.00 0.00 0.04 0.00 0.00 0.00 0.62
#&gt; SP101634     8  0.4858     0.0656 0.00 0.00 0.38 0.00 0.00 0.22 0.00 0.40 0.00
#&gt; SP101642     9  0.3408     0.6463 0.00 0.16 0.06 0.00 0.02 0.00 0.00 0.00 0.76
#&gt; SP101648     3  0.2537     0.6221 0.00 0.00 0.84 0.02 0.00 0.12 0.00 0.02 0.00
#&gt; SP101654     9  0.3637     0.6301 0.00 0.24 0.02 0.00 0.04 0.00 0.00 0.00 0.70
#&gt; SP101658     6  0.4105     0.5686 0.12 0.00 0.00 0.02 0.00 0.72 0.12 0.02 0.00
#&gt; SP101662     8  0.6791     0.1861 0.00 0.00 0.14 0.22 0.00 0.26 0.06 0.32 0.00
#&gt; SP101666     6  0.4635     0.1610 0.00 0.00 0.20 0.00 0.00 0.52 0.00 0.28 0.00
#&gt; SP101670     9  0.3347     0.4584 0.00 0.02 0.00 0.00 0.02 0.00 0.24 0.00 0.72
#&gt; SP101674     2  0.3919     0.0190 0.00 0.56 0.00 0.00 0.06 0.00 0.00 0.00 0.38
#&gt; SP101678     9  0.3402     0.6073 0.00 0.28 0.00 0.00 0.04 0.00 0.00 0.00 0.68
#&gt; SP101682     2  0.1473     0.6155 0.00 0.92 0.00 0.00 0.02 0.00 0.00 0.00 0.06
#&gt; SP101686     9  0.3014     0.6037 0.00 0.04 0.00 0.00 0.08 0.00 0.06 0.00 0.82
#&gt; SP101690     9  0.3789     0.3468 0.00 0.46 0.00 0.00 0.04 0.00 0.00 0.00 0.50
#&gt; SP101694     2  0.3572     0.5221 0.00 0.70 0.00 0.00 0.22 0.00 0.00 0.00 0.08
#&gt; SP101700     6  0.5772     0.1910 0.00 0.00 0.20 0.00 0.00 0.52 0.18 0.04 0.06
#&gt; SP101708     3  0.2537     0.6221 0.00 0.00 0.84 0.02 0.00 0.12 0.00 0.02 0.00
#&gt; SP101716     9  0.3712     0.4905 0.00 0.06 0.30 0.00 0.00 0.00 0.00 0.00 0.64
#&gt; SP101724     9  0.3402     0.6073 0.00 0.28 0.00 0.00 0.04 0.00 0.00 0.00 0.68
#&gt; SP101732     9  0.3316     0.6214 0.00 0.26 0.00 0.00 0.04 0.00 0.00 0.00 0.70
#&gt; SP101740     9  0.2821     0.6412 0.00 0.22 0.00 0.00 0.02 0.00 0.00 0.00 0.76
#&gt; SP101795     9  0.3221     0.6314 0.00 0.24 0.00 0.00 0.04 0.00 0.00 0.00 0.72
#&gt; SP101845     9  0.3116     0.6390 0.00 0.22 0.00 0.00 0.04 0.00 0.00 0.00 0.74
#&gt; SP101881     3  0.2537     0.6221 0.00 0.00 0.84 0.02 0.00 0.12 0.00 0.02 0.00
#&gt; SP101891     9  0.3221     0.6314 0.00 0.24 0.00 0.00 0.04 0.00 0.00 0.00 0.72
#&gt; SP101921     7  0.2890     0.7818 0.00 0.00 0.00 0.00 0.08 0.00 0.80 0.00 0.12
#&gt; SP101931     2  0.3193     0.2610 0.00 0.68 0.00 0.00 0.02 0.00 0.00 0.00 0.30
#&gt; SP102015     9  0.3316     0.6214 0.00 0.26 0.00 0.00 0.04 0.00 0.00 0.00 0.70
#&gt; SP102035     8  0.4886     0.1667 0.00 0.00 0.36 0.02 0.00 0.14 0.00 0.48 0.00
#&gt; SP102045     2  0.3478     0.4840 0.00 0.66 0.00 0.00 0.30 0.00 0.00 0.00 0.04
#&gt; SP102055     3  0.5023     0.2365 0.00 0.00 0.52 0.06 0.00 0.10 0.00 0.32 0.00
#&gt; SP102064     6  0.5620     0.1325 0.00 0.00 0.24 0.22 0.00 0.46 0.00 0.08 0.00
#&gt; SP102074     9  0.3221     0.6314 0.00 0.24 0.00 0.00 0.04 0.00 0.00 0.00 0.72
#&gt; SP102084     3  0.2001     0.6152 0.00 0.00 0.84 0.00 0.00 0.16 0.00 0.00 0.00
#&gt; SP102090     9  0.5156     0.4469 0.00 0.10 0.18 0.00 0.00 0.06 0.00 0.04 0.62
#&gt; SP102096     9  0.3997     0.2885 0.00 0.46 0.00 0.00 0.06 0.00 0.00 0.00 0.48
#&gt; SP102103     9  0.3985     0.6335 0.00 0.22 0.00 0.00 0.08 0.00 0.02 0.00 0.68
#&gt; SP102113     3  0.2537     0.6221 0.00 0.00 0.84 0.02 0.00 0.12 0.00 0.02 0.00
#&gt; SP102123     9  0.2431     0.6499 0.00 0.16 0.00 0.00 0.02 0.00 0.00 0.00 0.82
#&gt; SP102133     9  0.4099     0.5986 0.00 0.20 0.00 0.00 0.16 0.00 0.00 0.00 0.64
#&gt; SP102143     8  0.5060     0.5183 0.00 0.00 0.14 0.00 0.00 0.16 0.06 0.62 0.02
#&gt; SP102161     5  0.1670     0.7387 0.00 0.00 0.00 0.00 0.88 0.00 0.12 0.00 0.00
#&gt; SP102168     3  0.3221     0.5201 0.00 0.00 0.72 0.00 0.00 0.04 0.00 0.00 0.24
#&gt; SP102174     2  0.3545     0.1931 0.00 0.64 0.00 0.00 0.04 0.00 0.00 0.00 0.32
#&gt; SP102187     3  0.2275     0.6086 0.02 0.00 0.84 0.00 0.00 0.14 0.00 0.00 0.00
#&gt; SP102485     1  0.0000     0.9107 1.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; SP102489     7  0.1707     0.7421 0.08 0.00 0.00 0.02 0.00 0.00 0.90 0.00 0.00
#&gt; SP102499     6  0.4477     0.6230 0.26 0.00 0.00 0.18 0.00 0.56 0.00 0.00 0.00
#&gt; SP102507     8  0.2929     0.6612 0.00 0.00 0.00 0.24 0.00 0.02 0.00 0.74 0.00
#&gt; SP102511     6  0.4985     0.6321 0.22 0.00 0.00 0.18 0.04 0.56 0.00 0.00 0.00
#&gt; SP102517     5  0.1033     0.7890 0.06 0.00 0.00 0.00 0.94 0.00 0.00 0.00 0.00
#&gt; SP102523     7  0.1269     0.7869 0.00 0.00 0.00 0.00 0.08 0.00 0.92 0.00 0.00
#&gt; SP102529     5  0.1033     0.7890 0.06 0.00 0.00 0.00 0.94 0.00 0.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2654 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-8-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-8-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-8-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-9'>
<p><a id='tab-ATC-skmeans-get-classes-9-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 10), get_membership(res, k = 10))
</code></pre>

<pre><code>#&gt;          class entropy silhouette   p1   p2   p3   p4   p5   p6   p7   p8   p9  p10
#&gt; SP1003       2  0.2121    0.63721 0.00 0.88 0.00 0.00 0.02 0.00 0.00 0.00 0.06 0.04
#&gt; SP10084      2  0.0729    0.65688 0.00 0.96 0.00 0.00 0.00 0.00 0.00 0.00 0.04 0.00
#&gt; SP1009       8  0.0426    0.64740 0.00 0.00 0.00 0.02 0.00 0.00 0.00 0.98 0.00 0.00
#&gt; SP10150      2  0.5540    0.45244 0.00 0.54 0.00 0.00 0.18 0.00 0.10 0.00 0.14 0.04
#&gt; SP101515     3  0.4791    0.29702 0.00 0.00 0.50 0.28 0.00 0.00 0.02 0.20 0.00 0.00
#&gt; SP101519     7  0.4814    0.53397 0.00 0.04 0.00 0.00 0.04 0.00 0.60 0.00 0.24 0.08
#&gt; SP101521     9  0.3319    0.44752 0.00 0.30 0.00 0.00 0.00 0.00 0.00 0.00 0.66 0.04
#&gt; SP101523     3  0.0986    0.69031 0.00 0.00 0.94 0.00 0.00 0.06 0.00 0.00 0.00 0.00
#&gt; SP101526     2  0.2047    0.59509 0.00 0.82 0.00 0.00 0.00 0.00 0.00 0.00 0.18 0.00
#&gt; SP101528     8  0.5620    0.31684 0.00 0.00 0.10 0.34 0.00 0.00 0.02 0.42 0.00 0.12
#&gt; SP101532     8  0.6574    0.31277 0.00 0.02 0.10 0.12 0.02 0.24 0.00 0.44 0.00 0.06
#&gt; SP101536     5  0.4702    0.22964 0.00 0.16 0.00 0.00 0.50 0.00 0.00 0.00 0.32 0.02
#&gt; SP101540     2  0.3958    0.15896 0.00 0.52 0.00 0.00 0.02 0.00 0.00 0.00 0.42 0.04
#&gt; SP101544     2  0.1759    0.62384 0.00 0.86 0.00 0.00 0.00 0.00 0.00 0.00 0.14 0.00
#&gt; SP101548     8  0.4511    0.50185 0.00 0.00 0.06 0.28 0.00 0.04 0.02 0.60 0.00 0.00
#&gt; SP101552     9  0.3638    0.40610 0.00 0.28 0.00 0.00 0.02 0.00 0.00 0.00 0.66 0.04
#&gt; SP101558     3  0.1152    0.69336 0.02 0.00 0.94 0.00 0.00 0.04 0.00 0.00 0.00 0.00
#&gt; SP101564     9  0.2972    0.43752 0.00 0.28 0.00 0.00 0.00 0.00 0.00 0.00 0.70 0.02
#&gt; SP101572     2  0.2535    0.64271 0.00 0.84 0.00 0.00 0.10 0.00 0.00 0.00 0.02 0.04
#&gt; SP101576     9  0.3114    0.36614 0.00 0.32 0.00 0.00 0.00 0.00 0.00 0.00 0.66 0.02
#&gt; SP101580     9  0.2575    0.48690 0.00 0.28 0.00 0.00 0.00 0.00 0.00 0.00 0.72 0.00
#&gt; SP101584     2  0.3155    0.58994 0.00 0.78 0.00 0.00 0.00 0.00 0.04 0.00 0.14 0.04
#&gt; SP101588     9  0.4314    0.07887 0.00 0.02 0.00 0.00 0.00 0.00 0.40 0.00 0.50 0.08
#&gt; SP101592     2  0.3488    0.37063 0.00 0.60 0.00 0.00 0.00 0.00 0.00 0.00 0.36 0.04
#&gt; SP101596     2  0.3265    0.56534 0.00 0.74 0.00 0.00 0.02 0.00 0.00 0.00 0.20 0.04
#&gt; SP101600     3  0.1152    0.69336 0.02 0.00 0.94 0.00 0.00 0.04 0.00 0.00 0.00 0.00
#&gt; SP101604     8  0.5361    0.33749 0.00 0.00 0.16 0.38 0.00 0.04 0.02 0.40 0.00 0.00
#&gt; SP101610     4  0.6744   -0.28403 0.00 0.02 0.26 0.36 0.02 0.04 0.02 0.24 0.00 0.04
#&gt; SP101616     9  0.3165    0.53338 0.00 0.00 0.26 0.00 0.00 0.00 0.00 0.00 0.70 0.04
#&gt; SP101622     9  0.2579    0.54974 0.00 0.20 0.00 0.00 0.00 0.00 0.00 0.00 0.78 0.02
#&gt; SP101628     9  0.2692    0.53393 0.00 0.22 0.00 0.00 0.00 0.00 0.00 0.00 0.76 0.02
#&gt; SP101634     8  0.6880   -0.03969 0.00 0.02 0.32 0.12 0.02 0.14 0.00 0.32 0.00 0.06
#&gt; SP101642     9  0.2340    0.62999 0.00 0.02 0.08 0.00 0.00 0.00 0.00 0.00 0.86 0.04
#&gt; SP101648     3  0.0729    0.69462 0.00 0.00 0.96 0.00 0.00 0.00 0.00 0.04 0.00 0.00
#&gt; SP101654     9  0.2122    0.62869 0.00 0.10 0.04 0.00 0.00 0.00 0.00 0.00 0.86 0.00
#&gt; SP101658     6  0.4980    0.52530 0.06 0.00 0.08 0.12 0.00 0.66 0.06 0.00 0.00 0.02
#&gt; SP101662     8  0.5429    0.32614 0.00 0.00 0.18 0.36 0.00 0.04 0.02 0.40 0.00 0.00
#&gt; SP101666     8  0.6927   -0.02398 0.00 0.02 0.28 0.12 0.02 0.24 0.00 0.28 0.00 0.04
#&gt; SP101670     9  0.4679    0.04941 0.00 0.02 0.00 0.00 0.02 0.00 0.40 0.00 0.48 0.08
#&gt; SP101674     9  0.3973    0.00642 0.00 0.44 0.00 0.00 0.02 0.00 0.00 0.00 0.50 0.04
#&gt; SP101678     9  0.1594    0.61978 0.00 0.12 0.00 0.00 0.00 0.00 0.00 0.00 0.88 0.00
#&gt; SP101682     2  0.1412    0.64486 0.00 0.90 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.00
#&gt; SP101686     9  0.2870    0.57867 0.00 0.04 0.00 0.00 0.02 0.00 0.14 0.00 0.80 0.00
#&gt; SP101690     9  0.3224    0.29819 0.00 0.36 0.00 0.00 0.02 0.00 0.00 0.00 0.62 0.00
#&gt; SP101694     2  0.4748    0.38952 0.00 0.54 0.00 0.00 0.14 0.00 0.00 0.00 0.28 0.04
#&gt; SP101700     6  0.7661   -0.03569 0.00 0.02 0.26 0.12 0.02 0.26 0.20 0.08 0.00 0.04
#&gt; SP101708     3  0.0729    0.69462 0.00 0.00 0.96 0.00 0.00 0.00 0.00 0.04 0.00 0.00
#&gt; SP101716     9  0.3559    0.51774 0.00 0.00 0.26 0.00 0.00 0.00 0.02 0.00 0.68 0.04
#&gt; SP101724     9  0.1594    0.61978 0.00 0.12 0.00 0.00 0.00 0.00 0.00 0.00 0.88 0.00
#&gt; SP101732     9  0.1594    0.61978 0.00 0.12 0.00 0.00 0.00 0.00 0.00 0.00 0.88 0.00
#&gt; SP101740     9  0.1211    0.63484 0.00 0.08 0.00 0.00 0.00 0.00 0.00 0.00 0.92 0.00
#&gt; SP101795     9  0.1211    0.63484 0.00 0.08 0.00 0.00 0.00 0.00 0.00 0.00 0.92 0.00
#&gt; SP101845     9  0.1412    0.62986 0.00 0.10 0.00 0.00 0.00 0.00 0.00 0.00 0.90 0.00
#&gt; SP101881     3  0.0729    0.69462 0.00 0.00 0.96 0.00 0.00 0.00 0.00 0.04 0.00 0.00
#&gt; SP101891     9  0.1412    0.62730 0.00 0.10 0.00 0.00 0.00 0.00 0.00 0.00 0.90 0.00
#&gt; SP101921     7  0.3867    0.66549 0.00 0.02 0.00 0.00 0.10 0.00 0.74 0.00 0.10 0.04
#&gt; SP101931     2  0.2954    0.18467 0.00 0.58 0.00 0.00 0.00 0.00 0.00 0.00 0.42 0.00
#&gt; SP102015     9  0.1412    0.62986 0.00 0.10 0.00 0.00 0.00 0.00 0.00 0.00 0.90 0.00
#&gt; SP102035     8  0.6437   -0.06147 0.00 0.02 0.36 0.12 0.02 0.04 0.02 0.38 0.00 0.04
#&gt; SP102045     2  0.4921    0.40576 0.00 0.52 0.00 0.00 0.24 0.00 0.00 0.00 0.20 0.04
#&gt; SP102055     3  0.4820    0.27557 0.00 0.00 0.50 0.24 0.00 0.02 0.00 0.24 0.00 0.00
#&gt; SP102064     3  0.5240    0.31088 0.00 0.00 0.46 0.32 0.00 0.12 0.00 0.10 0.00 0.00
#&gt; SP102074     9  0.1412    0.62986 0.00 0.10 0.00 0.00 0.00 0.00 0.00 0.00 0.90 0.00
#&gt; SP102084     3  0.0850    0.69670 0.00 0.00 0.96 0.02 0.00 0.02 0.00 0.00 0.00 0.00
#&gt; SP102090     9  0.4297    0.44010 0.00 0.02 0.12 0.00 0.00 0.00 0.00 0.00 0.60 0.26
#&gt; SP102096     9  0.3649    0.21965 0.00 0.38 0.00 0.00 0.02 0.00 0.00 0.00 0.58 0.02
#&gt; SP102103     9  0.2586    0.62732 0.00 0.08 0.00 0.00 0.06 0.00 0.02 0.00 0.84 0.00
#&gt; SP102113     3  0.0729    0.69462 0.00 0.00 0.96 0.00 0.00 0.00 0.00 0.04 0.00 0.00
#&gt; SP102123     9  0.1955    0.62041 0.00 0.06 0.00 0.00 0.00 0.00 0.00 0.00 0.88 0.06
#&gt; SP102133     9  0.2947    0.59075 0.00 0.10 0.00 0.00 0.12 0.00 0.00 0.00 0.78 0.00
#&gt; SP102143     8  0.6627    0.35145 0.00 0.02 0.14 0.12 0.02 0.06 0.08 0.52 0.00 0.04
#&gt; SP102161     5  0.0986    0.76860 0.00 0.00 0.00 0.00 0.94 0.00 0.06 0.00 0.00 0.00
#&gt; SP102168     3  0.2489    0.50045 0.00 0.00 0.74 0.00 0.00 0.00 0.00 0.00 0.26 0.00
#&gt; SP102174     9  0.2996    0.07534 0.00 0.46 0.00 0.00 0.00 0.00 0.00 0.00 0.54 0.00
#&gt; SP102187     3  0.1152    0.69336 0.02 0.00 0.94 0.00 0.00 0.04 0.00 0.00 0.00 0.00
#&gt; SP102485     1  0.0000    0.88326 1.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt;  [ reached &#39;max&#39; / getOption(&quot;max.print&quot;) -- omitted 2661 rows ]
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-9-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-9-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-9-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-ATC-skmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-9'>k = 10</a></li>
</ul>
<div id='tab-ATC-skmeans-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-1"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-2"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-3"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-4'>
<pre><code class="r">consensus_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-4"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-5'>
<pre><code class="r">consensus_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-5"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-6'>
<pre><code class="r">consensus_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-6-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-6"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-7'>
<pre><code class="r">consensus_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-7-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-7"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-8'>
<pre><code class="r">consensus_heatmap(res, k = 9)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-8-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-8"/></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-9'>
<pre><code class="r">consensus_heatmap(res, k = 10)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-9-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-9"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-ATC-skmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-membership-heatmap'>
<ul>
<li><a href='#tab-ATC-skmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-9'>k = 10</a></li>
</ul>
<div id='tab-ATC-skmeans-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-1"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-2"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-3"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-4'>
<pre><code class="r">membership_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-4"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-5'>
<pre><code class="r">membership_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-5"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-6'>
<pre><code class="r">membership_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-6-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-6"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-7'>
<pre><code class="r">membership_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-7-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-7"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-8'>
<pre><code class="r">membership_heatmap(res, k = 9)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-8-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-8"/></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-9'>
<pre><code class="r">membership_heatmap(res, k = 10)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-9-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-9"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-signatures'>
<ul>
<li><a href='#tab-ATC-skmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-9'>k = 10</a></li>
</ul>
<div id='tab-ATC-skmeans-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-1-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-1"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-2-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-2"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-3-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-3"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-4'>
<pre><code class="r">get_signatures(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-4-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-4"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-5'>
<pre><code class="r">get_signatures(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-5-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-5"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-6'>
<pre><code class="r">get_signatures(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-6-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-6"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-7'>
<pre><code class="r">get_signatures(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-7-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-7"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-8'>
<pre><code class="r">get_signatures(res, k = 9)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-8-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-8"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-9'>
<pre><code class="r">get_signatures(res, k = 10)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-9-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-9"/></p>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-9'>k = 10</a></li>
</ul>
<div id='tab-ATC-skmeans-get-signatures-no-scale-1'>
<pre><code class="r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-1-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-1"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-2'>
<pre><code class="r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-2-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-2"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-3'>
<pre><code class="r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-3-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-3"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-4'>
<pre><code class="r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-4-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-4"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-5'>
<pre><code class="r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-5-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-5"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-6'>
<pre><code class="r">get_signatures(res, k = 7, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-6-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-6"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-7'>
<pre><code class="r">get_signatures(res, k = 8, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-7-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-7"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-8'>
<pre><code class="r">get_signatures(res, k = 9, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-8-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-8"/></p>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-9'>
<pre><code class="r">get_signatures(res, k = 10, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-ATC-skmeans-get-signatures-no-scale-9-1.png" alt="plot of chunk tab-ATC-skmeans-get-signatures-no-scale-9"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

```
#> Error in fit_diagram(combinations, "euler", input, shape, control, ...): !any(duplicated(names(combinations))) is not TRUE
```

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```




t-SNE plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-ATC-skmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-dimension-reduction'>
<ul>
<li><a href='#tab-ATC-skmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-5'>k = 6</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-6'>k = 7</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-7'>k = 8</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-8'>k = 9</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-9'>k = 10</a></li>
</ul>
<div id='tab-ATC-skmeans-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-4'>
<pre><code class="r">dimension_reduction(res, k = 5, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-5'>
<pre><code class="r">dimension_reduction(res, k = 6, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-6'>
<pre><code class="r">dimension_reduction(res, k = 7, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-7'>
<pre><code class="r">dimension_reduction(res, k = 8, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-8'>
<pre><code class="r">dimension_reduction(res, k = 9, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-9'>
<pre><code class="r">dimension_reduction(res, k = 10, method = &quot;t-SNE&quot;)
</code></pre>

<pre><code>#&gt; Error in Rtsne.default(X = structure(c(-3.3804321739366, -4.96756462224444, : Remove duplicates before running TSNE.
</code></pre>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk ATC-skmeans-collect-classes](figure_cola/ATC-skmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/functional_enrichment.html) for more detailed explanations.


 

## Session info


```r
sessionInfo()
```

```
#> R version 4.0.2 (2020-06-22)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS High Sierra 10.13.6
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8
#> 
#> attached base packages:
#> [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] genefilter_1.70.0         ComplexHeatmap_2.7.2.1000 markdown_1.1             
#>  [4] knitr_1.30                forcats_0.5.0             stringr_1.4.0            
#>  [7] dplyr_1.0.2               purrr_0.3.4               readr_1.4.0              
#> [10] tidyr_1.1.2               tibble_3.0.4              ggplot2_3.3.3            
#> [13] tidyverse_1.3.0           cola_1.9.0.1013           synchronicity_1.3.5      
#> [16] bigmemory_4.5.36          Biobase_2.48.0            BiocGenerics_0.34.0      
#> [19] roxytest_0.0.1            pacman_0.5.1             
#> 
#> loaded via a namespace (and not attached):
#>   [1] readxl_1.3.1         uuid_0.1-4           backports_1.2.1      circlize_0.4.11     
#>   [5] NMF_0.23.0           plyr_1.8.6           polylabelr_0.2.0     splines_4.0.2       
#>   [9] listenv_0.8.0        usethis_2.0.0        gridBase_0.4-7       digest_0.6.27       
#>  [13] foreach_1.5.1        htmltools_0.5.0      magick_2.5.2         fansi_0.4.1         
#>  [17] magrittr_2.0.1       memoise_1.1.0        cluster_2.1.0        doParallel_1.0.16   
#>  [21] openxlsx_4.2.3       remotes_2.2.0        annotate_1.66.0      globals_0.14.0      
#>  [25] modelr_0.1.8         matrixStats_0.57.0   languageserver_0.3.9 prettyunits_1.1.1   
#>  [29] colorspace_2.0-0     blob_1.2.1           rvest_0.3.6          ggrepel_0.9.0       
#>  [33] haven_2.3.1          xfun_0.19            RCurl_1.98-1.2       jsonlite_1.7.2      
#>  [37] callr_3.5.1          crayon_1.3.4         microbenchmark_1.4-7 bigmemory.sri_0.1.3 
#>  [41] roxygen2_7.1.1       impute_1.62.0        survival_3.2-7       brew_1.0-6          
#>  [45] iterators_1.0.13     glue_1.4.2           polyclip_1.10-0      registry_0.5-1      
#>  [49] gtable_0.3.0         GetoptLong_1.0.5     car_3.0-10           pkgbuild_1.2.0      
#>  [53] shape_1.4.5          abind_1.4-5          scales_1.1.1         DBI_1.1.0           
#>  [57] rngtools_1.5         rstatix_0.6.0        Rcpp_1.0.5           xtable_1.8-4        
#>  [61] clue_0.3-58          bit_4.0.4            foreign_0.8-81       mclust_5.4.7        
#>  [65] stats4_4.0.2         DT_0.17              htmlwidgets_1.5.3    httr_1.4.2          
#>  [69] RColorBrewer_1.1-2   ellipsis_0.3.1       factoextra_1.0.7     XML_3.99-0.5        
#>  [73] pkgconfig_2.0.3      dbplyr_2.0.0         AnnotationDbi_1.50.3 tidyselect_1.1.0    
#>  [77] rlang_0.4.10         reshape2_1.4.4       munsell_0.5.0        cellranger_1.1.0    
#>  [81] tools_4.0.2          cli_2.2.0            RSQLite_2.2.1        generics_0.1.0      
#>  [85] devtools_2.3.2       broom_0.7.3          evaluate_0.14        yaml_2.2.1          
#>  [89] bit64_4.0.5          processx_3.4.5       fs_1.5.0             zip_2.1.1           
#>  [93] future_1.21.0        slam_0.1-48          xml2_1.3.2           compiler_4.0.2      
#>  [97] rstudioapi_0.13      curl_4.3             png_0.1-7            testthat_3.0.1      
#> [101] ggsignif_0.6.0       reprex_0.3.0         stringi_1.5.3        highr_0.8           
#> [105] ps_1.5.0             desc_1.2.0           lattice_0.20-41      Matrix_1.3-2        
#> [109] vctrs_0.3.6          pillar_1.4.7         lifecycle_0.2.0      furrr_0.2.1         
#> [113] BiocManager_1.30.10  eulerr_6.1.0         GlobalOptions_0.1.2  bitops_1.0-6        
#> [117] data.table_1.13.6    R6_2.5.0             rio_0.5.16           IRanges_2.22.2      
#> [121] parallelly_1.23.0    sessioninfo_1.1.1    codetools_0.2-18     assertthat_0.2.1    
#> [125] pkgload_1.1.0        pkgmaker_0.32.2      rprojroot_2.0.2      rjson_0.2.20        
#> [129] withr_2.3.0          S4Vectors_0.26.1     hms_0.5.3            rmarkdown_2.6       
#> [133] rvcheck_0.1.8        carData_3.0-4        sigminer_1.2.0.99    Rtsne_0.15          
#> [137] skmeans_0.2-13       Cairo_1.5-12.2       ggpubr_0.4.0         lubridate_1.7.9.2
```



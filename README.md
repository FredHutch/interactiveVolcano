## Interactive Volcano Plot

An interactive shiny app for creating and editing volcano plots. A volcano plot is a type of scatter plot represents differential expression of features (genes for example): on the x-axis we typically find the fold change and on the y-axis the p-value.

### Run 

This shiny app takes a single file input. It is a dataframe with at least three columns: gene ID, significance value, and effect size. The dataframe shoudl be saved as either a `.csv`, `.tsv`, or `.txt`. 

To run the application, clone this repository and place your data file within the repository. Open app.R in RStudio and click Run App.

### Usage

Once running the application opens to a tab with your volcano plot. Hover over the plot points to view geneID and other metrics.

On the left hand sidebar you'll find various ways to cuostmize and annotate your plot including setting the axes variables, coloring the plot by differentially expressed gene, and labeling specific genes.

Easily download your volcano plot as a `.pdf` by clicking the download button.

![](/assets/volcanoPlotScreenShot.png)

On the second tab you'll find a rendered data table of the uploaded dataset that can be filtered to only show differentially expressed genes.

![](/assets/dataScreenShot.png)

## Interactive Volcano Plot

An interactive shiny app for creating and editing volcano plots. A volcano plot is a type of scatter plot represents differential expression of features (genes for example): on the x-axis we typically find the fold change and on the y-axis the p-value.

### Set up

The easiest way to install this application is to clone it from this GitHub. Open the command line (terminal on Mac) and type the following commands:

Go to the directory that you want to download the app to.
```
cd <PATH/TO/DIR>
```

Clone this repository
```
git clone https://github.com/FredHutch/interactiveVolcano.git
```

This app is built to work in a specific private Docker environment. To run locally, open `app.R` and switch the variable `local` to `local <- TRUE` 

### Adding your data 

This shiny app takes a single file input. It is a dataframe with at least three columns: gene ID, significance value, and effect size. The dataframe should be saved as either a `.csv`, `.tsv`, or `.txt`.

Replace the example dataset in `interactiveVolcano/data` with your dataset.

### Running the app

Run the app from RStudio:

1. Open `app.R`.
2. Click `Run App` in the top right hand corner.

>Make sure that RStudio's working directory points to the application

Run from the command line:

1. Navigate to the `interactiveVolcano` directory
2. `R -e "shiny::runApp('./')"`

### Application features

Once running the application opens to a tab with your volcano plot. Hover over the plot points to view geneID and other metrics.

On the left hand sidebar you'll find various ways to cuostmize and annotate your plot including setting the axes variables, coloring the plot by differentially expressed gene, and labeling specific genes.

Easily download your volcano plot as a `.pdf` by clicking the download button.

![](/assets/volcanoPlotScreenShot.png)

On the second tab you'll find a rendered data table of the uploaded dataset that can be filtered to only show differentially expressed genes.

![](/assets/dataScreenShot.png)

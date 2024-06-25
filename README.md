# SURE_example_3
 
In this tutorial, we showcase the application of SURE in tissues other than the hematopoietic system, such as the brain. We used the neuron data from the Human Brain Atlas. Users can download this data from the [Chan Zuckerberg CELL by GENE Discover](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443) or use the wget command provided below:

```bash
wget -c https://datasets.cellxgene.cziscience.com/b9171f05-8112-4a55-95f2-4cf8a57df8a2.h5ad
```

Here are the code scripts for the analysis. Users can adjust the code according to their running environment.
1. [Prepare the input files for computing primary metacells.](./batch_prepare_mtx_files_4_SURE.py) Since the neuron compartment of the Human Brain Atlas contains nearly 2.5 million cells, this code will divide the data into several batches, with each batch containing approximately 10,000 cells.
2. [Compute the primary metacells.](./HBCA_batch_SURE_train.py) This is a batch processing file that automatically processes the data batches obtained through splitting in the previous step.
3. [Aggregate the expression profiles of cells to calculate the expression profiles of primary metacells.](./batch_prepare_primary_metacells.py)
4. [Prepare the input files for computing secondary metacells.](./prepare_secondary_mtx.ipynb)
5. [Compute the secondary metacells.](./HBCA_SURE_secondary_train.ipynb)
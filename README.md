# LoonDataCollab
Loon data exploration and predictive networking

The datasets are too large to upload to github (even after compressing). 

You will need three files: 
1) Download the file 'loon-flights-2021Q1' from: https://zenodo.org/record/5119968#.ZFqPPoTMLZQ for working with GPS data
2) Download the file 'link_reports-20210101-00' and 'link_intents' from: 'https://zenodo.org/record/6629754' for working with network data

### Notebooks: 
- Loon Aerospace Data Exploration v2: uses GPS data to compute distances between a set of nodes 
- Link Capacities: uses network data to compute data rates between a set of nodes 
- Loon Data Genie: uses the link intents file to evaluate which links between nodes were made by the TS-SDN 
- DataVis: visualize nodes, links, and distances. Needs map.json
- Topology Selection Baseline: Explores topology selection problem using or-tools

### Python Files: 
- CommonImports.py: used to import common python libraries
- TopologySelectionAutomation.py: Contains functions used in topology selection baseline notebook

### Packages: 
- numpy v1.21.5
- pandas v1.4.2
- scipy v1.10.1
- matplotlib v3.5.1
- KeplerGl v0.3.2
- or-tools: v9.6.2543
- json

## Notes on DataVis and KeplerGL
- KeplerGL user guide: https://docs.kepler.gl/docs/user-guides 
- KeplerGL Jupyter Notebook Documentation: https://docs.kepler.gl/docs/keplergl-jupyter 
- KeplerGL does not display data for NaN values, which might result in loons disappearing in the visualization
- It is currently only showing the links requested at that specific time, not the established links
- To make changes to the configuration file, load up the map on Jupyter Notebook, make the changes within the Jupyter Notebook and then run the code block below it to save the configuration file and then apply it to the HTML file
- The KeplerGL map/widget is often glitchy when displayed on the Jupyter Notebook, so you might need to re-run all the code for it to work 
-All the data needs to be within one data frame in KeplerGL for the filters (time playback feature) to work effectively

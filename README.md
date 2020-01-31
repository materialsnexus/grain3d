This is a 3D grain growth code which simulates normal grain growth.

The code was created by Dr Menachem Lazar with contributions from Jeremy Mason, David Srolovitz and Robert MacPherson. The code was edited and refined by Dr Jonathan Bean for public use. 

The code is based on the MacPherson-Srolovitz evolution of grains during grain growth described in the paper below. 

Lazar, E. A., Mason, J. K., MacPherson, R. D., & Srolovitz, D. J. (2011). A more accurate three-dimensional grain growth algorithm. Acta Materialia, 59(17), 6837–6847. https://doi.org/10.1016/j.actamat.2011.07.052

To install the code you can do the following:

    git clone 3d-grain-growth
    cd 3d-grain-growth
    make

This will make a binary (executable) file which can be run on windows or linux(unix). You make need to edit the make file to specify which compiler is being used. 

The code requies some initial vertices in order to evolve. Vertices are simple a 3 vector which describes the points which reprent the grain structure. There are however some test files in the root directory called "Sample1.dat", "Sample2.dat" and "Sample3.dat". 

You can run the code in the following way:

    ./3dgrain Sample1.dat grain3d.param

There are two input files required by the code. Ther first is a file contians some preamble and the vertices used for the evolution. The format is as follows

    Line1:				step_count
    Line2:				total_time
    Line3:				clocks
    Line4:				triangles_deleted
    Line5:				digons_deleted
    Line6:				edges_deleted
    Line7:				tetrahedra_deleted
    Line8:				footballs_deleted
    Line9:				edge_nodes_deleted
    Line10:				edge_nodes_added
    Line11:			    node_count
	Lines12-onwards:	The 3 vectors of the vertices
	Lines(11+node_count) After the vertices have been defined the connectivity of each grain should be defined

Most of the parameters specified in lines1-11 can be set to 0 unless restarting a calculation.

The .param file contains variables for how the system will evolve over time.

    Line1: end_body							(The simulation will continue until this many grains are remaining)
	Line2: refine_edges_number				(How many steps to wait before refining the edges)
	Line3: remove_system_edge_nodes_number  (How many steps to wait before removing systeme edge nodes)
	Line4: output_summary_number			(How many steps to wait before outputting summary statistics)
	Line5: output_whole_system_number		(How many steps to wait before outputting the whole system as a single file)
	Line6: output_specific_neigh			(How many steps to wait before outputting the whole system grain and its neighbours grain by grain)
	Line7: output_specific					(How many steps to wait before outputing the whole system grain by grain)

The distributed grain3d.param file gives typical values for all of the input parameters.

The output files are of the .fe format and can be read by the Surface Evolver program created by Ken Brakke.

Note that in the model it is assumed that all grain boundary types have the same grain boundary energy. In a real system this is likely to be true but no simple methodology for determining the distribution of grain boundary energy exists. Furthermore experimental evidence of the grain boundary character distribution is disputed.



If any publications arise using the code please cite the following paper.

Lazar, E. A., Mason, J. K., MacPherson, R. D., & Srolovitz, D. J. (2011). A more accurate three-dimensional grain growth algorithm. Acta Materialia, 59(17), 6837–6847. https://doi.org/10.1016/j.actamat.2011.07.052

Please do get in touch if you want to understand more, can't make it work, find any bugs, if you are interested in making significant changes.

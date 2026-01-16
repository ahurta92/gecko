Here I am collecting all the code that I want merge into Gecko.

All of the code that I want to migrate will be locatioed in the migration folder at the top level of the directory.


## parsing tools

Within migration we have collected several parsing tools.  Ultimately, we only want a single dalton.py and madness.py for parsing Dalton output or madness output json.  In daltonToJson, we have dalton parsers that we need to incorporate into dalton.py.  dalton.py only converts data by parsing from output and converting into the datatype of the property.


## typical workflow

We can see a typical workflow in make_csv_data.ipynb.  In this notebook we see tools to create hyperpolarizability databases from MADNESS and Dalton files and then creates some plots for analysis.  I want to standardize this procedure.

## migration.cli

in this directory we have some scripts for submitting dalton and madness raman calculations along with slurm scripts.  For this repo, we don't want to deal so much with database generation. But at least, we want to proivde general, and easy to use flexible tools for generating differnt types of dalton calculations.

## migration.db

In this directory we have developed several tools for creating madness and dalton calculations.  This was built for raman calculations, so specifically it was designed to take Dalton template inputs, with molecules, write the molecule dalton input, and then for a raman calculation, first submit the optimization calculation, parese the new geometry, and write the raman calculation with the new geometry.  We want to both simplify and generalize the idea.  The basic idea is that some properties (like raman) are formed from multiple part calculations, in which in step requires, one or several outputs from previous calculations.  Other calculations, like a polarizabilty calculation at a defined geometry would be a single point.  Like, I said, I don't want to hard stuck how to structure a geometry, but I would like to development "calculations" which can be a single step or multistep and of course use the parses to parse the data along the way.  The thing's that are missing here are the hyperpolarizability.  For the Dalton implementation, I think the easiest way I can see this working is by somehow as input, taking a list of template inputs, where the number of elements in the list define the number of steps and the the actual inputs will be the template dalton input files.

## migration.viz

Currently contains loots of tools for geometric analysis of hyperpolarizability tensor component basis set error using the unit-sphere representation.  Application.py defines a trame application which reads in shg_ijk.csv data and allows a user to inspect the unit-sphere representation in both mra reference and basis sets and allows to user to inspect useful error metrics, compared across basis sets.  Field error defines those metrics.  These are very useful tools.

## migration.ideas

Something that isn't currently implemented but I think would be very useful is to have a single source of truth in terms of different properties.  For example, for hyperpolarizability analysis, we need primarily a single csv or dataframe containing the molecule, basis, and the components at each input frequency.  What would be useful is a script which can use the output parses to read new calculation and add to the shg if the calculation is not available.

## ideas

We don't want to databases to be hard stuck in where calculations are located.  If we are working with legacy databases, all we want to do is provide some directory iterator functions and a list of the molecule's, basis set calculations that we are looking for and call the the parser on those output files.

There should be a read only options, which only looks for calculations  using what I am describing above and and then read/write version which basically computes a database with some predefined structure which enables look up.

Something that would be nice for my users is a GUI option where a user can select a series of calculations to run the parser on and then the parser will be able to provide dataframe version of the data for easy read in and analysis.

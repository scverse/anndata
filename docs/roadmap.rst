Roadmap
=======

Documentation
-------------

-  A more in-depth prose description of how objects can be combined
-  Update documentaion on capabilities of backed mode

Concatenation
-------------

-  Provide a general api for selecting which elements to combine, and
   how they are combined

On disk
-------

-  Further formalize the on disk schema
-  Unify the code for different on disk formats (e.g. ``hdf5`` and
   ``zarr``)
-  Add an interface for extending. Allow third party packages to
   register how their data-type should be written and read from disk.
-  Allow arbitrary python objects to be stored, via pickling them to
   binary blobs

Distributed / out of core
-------------------------

-  Figure out the right model for out-of-core workflows

   -  One option is to “roll our own” similar to ``loompy``
   -  Another would be to integrate with an existing tool like ``dask``

-  Add an API for partial reading of ``AnnData`` objects. Specify only a
   subset of the elements and/ or a set of variables/ observations to
   read in at once.

More array types
----------------

-  Make AnnData objects more generic about what kinds of arrays can be
   used. For example, it would be good if we could work with arrays
   from: ``cupy``, ``sparse``, ``xarray``, ``arrow``

Benchmarks
----------

-  Add benchmarks suite (likely using `asv <https://asv.readthedocs.io>`__) to the repository

   -  Figure out how to run these as a CI service

Multimodal (more axes)
----------------------

-  Expand the model to handle multiple sets of variables

   -  Potentially take an xarray ``axes`` type approach

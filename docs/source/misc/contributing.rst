.. _contributing:

============
Contributing
============

We welcome contributions to the dorado package in the form of pull requests and issues made in the source `repository <https://github.com/passaH2O/dorado>`_.


Issues
------

If you are having any problems using dorado we suggest opening an `issue <https://github.com/passaH2O/dorado/issues>`_.
When you open an issue, please provide a clear description of what your scenario is, and what error message you are receiving.
If possible, please include a minimal working example of some code that breaks, and the error output you receive.

If there is some functionality you would like to see added to dorado you can also open an issue up to discuss that.
This can be code you plan to write and contribute, or it can be something you would like to have available but are not comfortable coding yourself.
Either way we are happy to help!


Pull Requests
-------------

If you have a feature that you would like to propose be integrated into dorado, then you should open a `pull request <https://github.com/passaH2O/dorado/pulls>`_.
To create a pull request, we recommend first `forking <https://guides.github.com/activities/forking/>`_ the repository, and then creating a separate branch to develop your feature (for reference see the `GitHub flow guide <https://guides.github.com/introduction/flow/>`_ and this `Git branching <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_ guide).
Then you can commit and develop your feature in your branch.
We ask that new features be accompanied by additional unit tests, to ensure that they operate as expected.
For unit testing, we use `pytest <https://docs.pytest.org/en/latest/>`_.
When you are satisfied with the code you have developed, you can open a pull request to the `"master" branch <https://github.com/passaH2O/dorado/tree/master>`_ of the project repository.
Please write a concise and descriptive title for the pull request, and provide a clear description of what the feature does and why you are proposing its addition to the project.
If your pull request does not pass our `continuous integration <https://docs.github.com/en/actions/building-and-testing-code-with-continuous-integration/about-continuous-integration>`_ checks, then we will not be able to merge your addition into the code.
We use `GitHub Actions <https://github.com/features/actions>`_ to ensure that dorado compiles and installs on MacOS, Windows, and Ubuntu on all supported versions of Python.
Before developing any new code, you are more than welcome to open an issue first to discuss your proposed addition.


Code Style
----------

To ensure consistency within the codebase, we follow some standard Python conventions for our formatting.
As much as possible we try to stick to the `PEP-8 <https://www.python.org/dev/peps/pep-0008/>`_ standard for our code.
For docstrings, we try to follow the `PEP-257 <https://www.python.org/dev/peps/pep-0257/>`_ standard.

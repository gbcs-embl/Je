# HowTo update the toolshed repository

After linting/testing your update using ``planemo test`` and ``planemo serve`` you can upload the changes to the toolshed repos.

For this you need to have the credentials or api_key of our ``gbcs-embl-heidelberg`` user.
To make it easier, put this information in your global ``~/.planemo.yml`` file. The format is described here: https://planemo.readthedocs.io/en/latest/configuration.html

The commands below work due to the ``.shed.yml`` file in this directory.

## Update testtoolshed

``planemo shed_update --shed_target testtoolshed``

Look at the result and then update the main toolshed.

## Update toolshed

``planemo shed_update``


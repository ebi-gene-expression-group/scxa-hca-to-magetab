# scxa-hca-to-magetab
## Instructions to install and test run locally:
```
cd \your\local\working\dir
mkdir logs
git clone git@github.com:ebi-gene-expression-group/scxa-hca-to-magetab.git
export PYTHONPATH=
virtualenv -p python3.7 vhca2mtab
pip3.7 install -r scxa-hca-to-magetab/requirements.txt
cd scxa-hca-to-magetab
python hca2mtab.py test .
```
## Test run output
- N.B. The test run restricts the number of bundles retrieved from HCA DCC to maximum 500 per technology)
- The experimental metadata in MAGETAB format for the imported HCA projects can be found in \your\local\working\dir
- The log file for your test run in \your\local\working\dir\hca2mtab.test.YYYY-MM-DD.log
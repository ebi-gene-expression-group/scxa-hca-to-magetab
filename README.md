# scxa-hca-to-magetab
## Instructions to install and test run locally:
cd \your\local\working\dir
mkdir logs
git clone git@github.com:ebi-gene-expression-group/scxa-hca-to-magetab.git
export PYTHONPATH=
virtualenv -p python3.7 vhca2mtab
pip3.7 install -r scxa-hca-to-magetab/requirements.txt
cd scxa-hca-to-magetab
python hca2mtab.py test .
## OUTPUT
- The experimental metadata in MAGETAB format for the imported HCA projects can be found in \your\local\working\dir
- The log file for your test run in \your\local\working\dir\hca2mtab.test.YYYY-MM-DD.log
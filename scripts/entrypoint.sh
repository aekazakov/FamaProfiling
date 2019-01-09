#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  mkdir -p /data/famaprofiling
  cd /data/famaprofiling
  pwd

  echo "downloading background database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_background_db.faa.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_background_db.faa.gz
  gunzip /data/famaprofiling/fama_background_db.faa.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/fama_background_db.faa --db /data/famaprofiling/fama_background_db
  rm /data/famaprofiling/fama_background_db.faa

  echo "downloading nitrogen cycle database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz
  tar xvf /data/famaprofiling/fama_nitrogen.tar.gz
  rm /data/famaprofiling/fama_nitrogen.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/fama_nitrogen_db.faa --db /data/famaprofiling/fama_nitrogen_db
  rm /data/famaprofiling/fama_nitrogen_db.faa

  echo "downloading universal markers database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_universal.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_universal.tar.gz
  tar xvf /data/famaprofiling/fama_universal.tar.gz
  rm /data/famaprofiling/fama_universal.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/fama_universal_db.faa --db /data/famaprofiling/fama_universal_db
  rm /data/famaprofiling/fama_universal_db.faa

  echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz
  tar xvf /data/famaprofiling/fama_taxonomy.tar.gz
  rm /data/famaprofiling/fama_taxonomy.tar.gz

  if [ -s "/data/famaprofiling/fama_background_db.dmnd" ] ; then
#  if [ -s "/data/famaprofiling/fama_background_db.dmnd" -a -s "/data/diamondsearchexampledb/fama_universal_db.dmnd"] ; then
    echo "DATA DOWNLOADED SUCCESSFULLY"
    touch /data/__READY__
  else
    echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi

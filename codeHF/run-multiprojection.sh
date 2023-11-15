#!/bin/bash

CONFIG=multiprojection_config.json

for dl in `seq 0.000 0.002 0.020` ; do
  echo "dl ${dl}"
  sed -i -E "s/dl0\.[0-9]{3}/dl${dl}/g" ${CONFIG}
  python3 multiprojection.py ${CONFIG}
done

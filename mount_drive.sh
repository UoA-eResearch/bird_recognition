#!/bin/bash

drive_name="ressci202000004-Rifleman/2019 BSMI/"
share="//files.auckland.ac.nz/research/${drive_name}"
mountpoint="data"
common_options="iocharset=utf8,workgroup=uoa,uid=${USER},nodev,nosuid"
options="username=${USER},${common_options}"

mkdir -p ${mountpoint}
sudo mount -t cifs "${share}" "${mountpoint}" -o "${options}"
if [ "$?" -gt "0" ]; then
  rmdir ${mountpoint}
fi

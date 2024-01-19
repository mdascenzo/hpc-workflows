#!/bin/bash
FILE_TO_UPLOAD=$1
BUCKET_NAME=$2
S3_KEY=$3

aws s3 cp "${FILE_TO_UPLOAD}" "s3://${BUCKET_NAME}/${S3_KEY}"

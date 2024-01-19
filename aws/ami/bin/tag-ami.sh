#!/bin/bash

# Input arguments
AMI_ID=$1

# Tag the AMI with its own ID
aws ec2 create-tags --resources $AMI_ID --tags Key=parallelcluster:image_id,Value=$AMI_ID

# Tag the AMI as available
aws ec2 create-tags --resources $AMI_ID --tags Key=parallelcluster:build_status,Value=available
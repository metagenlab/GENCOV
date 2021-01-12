# Contribution Guide
Frst of all thanks for thinking about contributing to this project.

[[_TOC_]]

## Introduction

We are always glad to recieve help, in any way, shape or form. 
This means filing bug reports, requesting features, providing feedback 
or adding to the code base.

## Suggesting a feature

Keeping the scope of the pipeline limited was the overall goal, but if there
is a glaring omission, feel free to shoot us an email.  

## Submitting Code 

It is great that you are willing to submitt your own code. 
You can email us to get development access for a feature_branch. 

The use of [f-strings](https://realpython.com/python-f-strings/) currently 
requires python 3.6+ as a minimum.  

Keep in mind that we have two protected branches. 
**Master** and **development_unstable** that are both running CI

### Branch: master
The master branch schould only ever contain running code
and is mostly used as our release branch 

### Branch: development_unstable
Contains code that is being integrated from feature branches. 
Can contain non running code, that is in the process of being fixed. 

Please make sure that your local feature_branches pass at least the basic tests
before making a pull request for dev_unstable.   

## Issuing a new release

* Please ensure that the latest commits **pass** the quicktests (also in the continous 
  integration pipelines!).
* Verify with more comprehensive tests, that the outputs of the pipeline
  still make sense. (Check for silently failing snakemake rules e.g. unexpected empty files)
* Inform the other active members of this project  <a href="mailto:BI-Support@rki.de,drechselo@rki.de,kmiecinskir@rki.de,r.w.kmiecinski@gmail.com?subject=GITLAB-PROJECT: CoV-Pipe - new release">Email Us</a>

### Instructions

Please exectute the *new_release.sh* script in the project root
    
    >$ ./new_release [patch,minor,major]

This script will only work in the **master** and **development_unstable** branch. 
Before running, please ensure that you are up to date on both. 
It wil automatically incremnt the semantic version and edit the __version__.py 
file for you.

You will recieve a prompt, asking if you are sure that you want to release the new
version. 
If you reply positively, it automatically tags and merges the latest 
commit (after amending it with the version file), into the master branch.
Additionally a new recipe will be added to the .conda recipe folder.

If the release process encounters an error, 
you can continue where you left of, after fixing the issue, by simple re-executing

    >$ ./new_release

You are still required to add proper release notes [here](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe/-/tags)  


#### Example 
Minor Release
    
    >$ git pull && git fetch origin
    ...
    >$ ./new_release.sh minor
    Current Version: v0.4.0
    Trying to create new minor release! (v0.5.0)
    You are in branch: development_unstable
    You are about to release Version 0.5.0
    This is a minor release.
     
    Ensure that you have run all required tests
    and cooridinated with the other members of the project
 
    Do you still want to continue? [Y/n]

Major Release
    
    >$ git pull && git fetch origin
    ...
    >$ ./new_release.sh ma
    Current Version: v0.4.0
    Trying to create new major release! (v1.0.0)
    You are in branch: development_unstable
    You are about to release Version 1.0.0
    This is a major release.
 
    Ensure that you have run all required tests
    and cooridinated with the other members of the project
 
    Do you still want to continue? [Y/n]

# Thank You!!

For taking the time to read through our contribution guide. We hope that you are
in no way discuraged from contributing to this project. 



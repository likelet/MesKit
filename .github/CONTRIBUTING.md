# nf-core/multiexseq: Contributing Guidelines

Hi there! Many thanks for taking an interest in improving multiexseq.

We try to manage the required tasks for multiexseq using GitHub issues, you probably came to this page when creating one. Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome! Contributions to the code are even more welcome ;)


## Contribution workflow
If you'd like to write some code for multiexseq, the standard workflow
is as follows:

1. Check that there isn't already an issue about your idea in the
   [multiexseq issues](https://github.com/likelet/multiexseq/issues) to avoid
   duplicating work.
    * If there isn't one already, please create one so that others know you're working on this
2. Fork the [multiexseq repository](https://github.com/multiexseq) to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [basic docs from GitHub](https://help.github.com/articles/fork-a-repo/) or even their [excellent interactive tutorial](https://try.github.io/).


## Tests
When you create a pull request with changes, [Travis CI](https://travis-ci.org/) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:


### Pipeline Tests
Each nf-core pipeline should be set up with a minimal set of test-data.
Travis CI then runs the pipeline on this data to ensure that it exists successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of Nextflow and also the minimum required version that is stated in the pipeline code.

## Getting help
For further information/help, please consult the [multiexseq documentation](https://github.com/multiexseq#documentation) and don't hesitate to get in touch on [me](zhaoqi@sysucc.org.cn)

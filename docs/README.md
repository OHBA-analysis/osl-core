### Setup instructions

For local usage, install gems 

    bundle install

May first need

    brew install libxml2 libxslt
    sudo gem install nokogiri -v 1.6.8.1 -- --use-system-libraries=true --with-xml2-include=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk/usr/include/libxml2

Then run the local server

    ./local_server.sh

### Adding pages

Still in progress, but the workflow will likely be

1. Modify Matlab file with desired content
2. Run `osl_publish('filename.m')`
3. `git add .`, `git commit`, `git push`

And the live site should be updated immediately
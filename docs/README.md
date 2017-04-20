### Setup instructions

When developing the tutorial pages, it can be helpful to view the pages locally. First make sure you have all the dependencies 

	gem install bundler
    brew install libxml2 libxslt
    sudo gem install nokogiri -v 1.6.8.1 -- --use-system-libraries=true --with-xml2-include=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk/usr/include/libxml2

Then install the GitHub jekyll dependencies

    bundle install

Finally, run the local server

    ./local_server.sh

### Adding pages

1. Modify Matlab file with desired content
2. Run `osl_publish('filename.m')`. View locally using the local server if required
3. `git add .`, `git commit`, `git push`

And the live site should be updated immediately
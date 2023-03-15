all: serve

check: build
	bundle exec htmlproofer ./_site/ --only-4xx --disable-external --ignore-empty-alt

serve:
	bundle exec jekyll serve -s . -d _site/

build:
	bundle exec jekyll build -s . -d _site/

push:
	git commit -am "`date`" && git push

install:
	gem install jekyll
	gem install bundler
	bundle config set path 'vendor/bundle'
	bundle install

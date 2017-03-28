---
layout: default
---

The OHBA Software Library (OSL) is created by the OHBA Analysis Group, OHBA, Oxford, UK.

OSL is a set of Matlab tools and scripts for running M/EEG analysis on CTF and Elekta Neuromag data, and is written mainly by members of the OHBA Analysis Group. It uses a combination of FSL, SPM and Fieldtrip.

### Get a copy of OSL

The Github repository is insufficient to run OSL, as it requires a bunch of supporting directories. These can be downloaded from the OSL wiki. To get access to this please email [mark.woolrich@ohba.ox.ac.uk](mailto:mark.woolrich@ohba.ox.ac.uk)

### Examples

{% for cat in site.category-list %}
<ul>
  {% for page in site.pages %}
    {% if page.resource == true %}
      {% for pc in page.categories %}
        {% if pc == cat %}
          <li><a href="{{ site.baseurl }}{{ page.url }}">{{ page.title }}</a></li>
        {% endif %}   <!-- cat-match-p -->
      {% endfor %}  <!-- page-category -->
    {% endif %}   <!-- resource-p -->
  {% endfor %}  <!-- page -->
</ul>
{% endfor %}  <!-- cat -->

<!-- ### Tables

| head1        | head two          | three |
|:-------------|:------------------|:------|
| ok           | good swedish fish | nice  |
| out of stock | good and plenty   | nice  |
| ok           | good `oreos`      | hmm   |
| ok           | good `zoute` drop | yumm  |

### Small image

![](https://assets-cdn.github.com/images/icons/emoji/octocat.png)

### Large image

![](https://guides.github.com/activities/hello-world/branching.png)

### Code

```matlab
x = f(x)
``` -->

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree: .
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% set shown_methods = methods | reject("==", "__init__") | reject("in", inherited_members) | list %}
   {% if shown_methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: .
      {% for item in shown_methods %}
      ~{{ name }}.{{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

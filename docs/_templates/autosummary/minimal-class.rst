{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

..  autoclass:: {{ objname }}

    {% block attributes %}
    {% if attributes %}
    .. rubric:: Attributes

    {% for item in attributes %}
    .. autoattribute:: {{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block methods %}
    {% set shown_methods = methods | reject("==", "__init__") | reject("in", inherited_members) | list %}
    {% if shown_methods %}
    .. rubric:: Methods

    {% for item in shown_methods %}
    .. automethod:: {{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}

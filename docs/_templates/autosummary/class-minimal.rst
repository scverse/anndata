{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

..  autoclass:: {{ objname }}

    {% block attributes %}
    {%- for item in attributes %}
    {%- if loop.first %}
    .. rubric:: Attributes

    .. autoattribute:: {{ item }}
    {% endif %}
    {%- endfor %}
    {% endblock %}

    {% block methods %}
    {%- for item in methods if item != "__init__" and item not in inherited_members %}
    {%- if loop.first %}
    .. rubric:: Methods
    {% endif %}
    .. automethod:: {{ item }}
    {%- endfor %}
    {% endblock %}

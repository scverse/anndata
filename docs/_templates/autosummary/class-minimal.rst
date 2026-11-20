{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

..  autoclass:: {{ objname }}
    :show-inheritance:

    {% block attributes %}
    {%- for item in attributes %}
    {%- if loop.first %}
    .. rubric:: Attributes
    {% endif %}
    .. auto{{ [fullname, item] | join(".") | member_type }}:: {{ item }}
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

from typing import Iterable

def make_single_column_table(heading: str, entries: Iterable[str]) -> str:
    table = f"""
        <table>
            <thead>
                <tr><th> {heading} </th></tr>
            </thead>
            <tbody>
            """
    for entry in entries:
        table += f"""
                <tr>
                    <td>{entry}</td>
                </tr>
                """
    table +=  """
            </tbody>
        </table>
    """  
    return table

def make_single_row_table(entries: Iterable[str]) -> str:
    table = f"""
        <table>
            <tbody>
                <tr>
            """
    for entry in entries:
        table += f"""<td style="vertical-align:top">{entry}</td>"""
    table +=  """
                </tr>
            </tbody>
        </table>
    """
    return table

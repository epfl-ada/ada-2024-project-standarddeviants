def save_plotly(
    fig, savename: str, rgb_bgcolor: tuple[int, int, int] = (34, 37, 41)
) -> None:
    """Save plotly figure to html and svg formats for website

    Args:
        fig (plotly figure): plotly figure
        savename (str): path like without extensions
        rgb_bgcolor (tuple[int, int, int], optional): The RGB colors of the background of the website. Defaults to (34, 37, 41).
    """

    fig.update_layout(
        plot_bgcolor=f"rgb({rgb_bgcolor[0]}, {rgb_bgcolor[1]}, {rgb_bgcolor[2]})",
        paper_bgcolor=f"rgb({rgb_bgcolor[0]}, {rgb_bgcolor[1]}, {rgb_bgcolor[2]})",
    )

    custom_style = f"""
    <style>
    body, html {{
        background-color: rgb({rgb_bgcolor[0]}, {rgb_bgcolor[1]}, {rgb_bgcolor[2]});
    }}
    </style>
    """

    # Save the HTML file and inject the custom style
    html_content = fig.to_html(full_html=True, include_plotlyjs="cdn")
    # Insert the style in the <head>
    html_content = html_content.replace("</head>", f"{custom_style}</head>")

    # Write to a file
    with open(f"{savename}.html", "w") as f:
        f.write(html_content)

    fig.write_image(f"{savename}.svg")

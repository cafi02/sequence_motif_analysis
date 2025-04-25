import pandas as pd
import re


def filter_e_value_subset(df, e_value):
    """
    Filters a DataFrame based on an E-value threshold.

    Parameters:
        df (pandas.DataFrame): DataFrame containing an "E-value" column.
        e_value (float): Threshold for filtering E-values.

    Returns:
        pandas.DataFrame: Filtered DataFrame with rows having E-value ≤ e_value.
    """
    # Convert the column to numeric, forcing errors to NaN (if any non-numeric values exist)
    #df["Full Sequence: E-value"] = pd.to_numeric(df["Full Sequence: E-value"], errors='coerce')
    df.loc[:, "Full Sequence: E-value"] = pd.to_numeric(df["Full Sequence: E-value"], errors='coerce')

    # Filter rows where the E-value is ≤ the threshold
    result_df = df[df["Full Sequence: E-value"] <= e_value]

    print(f"Extracted {len(result_df)} rows for E-value <= {e_value}")

    return result_df

def split_phmmer_keywords(df, keywords_dict):
    """
    Splits the input DataFrame into multiple DataFrames based on keyword groups.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'description of target' column.
        keywords_dict (dict): A dictionary where:
            - Keys define the names of the result DataFrames.
            - Values are lists of keywords to search for in the 'description of target' column.

    Returns:
        dict: A dictionary where keys correspond to the categories in keywords_dict
              and values are DataFrames containing rows that match any keyword in that category.
              Also includes a 'remaining' DataFrame with unmatched rows.
    """
    description = df['Description']
    result_dfs = {key: pd.DataFrame(columns=df.columns) for key in keywords_dict}  # Initialize empty DataFrames

    matched_indices = set()  # Track matched row indices to avoid modifying df during iteration

    for category, keyword_list in keywords_dict.items():
        keyword_pattern = '|'.join(map(re.escape, keyword_list))  # Create regex pattern from keyword list
        matched_rows = df[description.str.contains(keyword_pattern, case=False, na=False, regex=True)]

        result_dfs[category] = matched_rows  # Store matched rows
        matched_indices.update(matched_rows.index)  # Keep track of matched indices

        print(f"Extracted {len(matched_rows)} rows for category: {category}")

    # Create DataFrame for remaining unmatched rows
    result_dfs['Remains'] = df.drop(matched_indices)

    print(f"Remaining rows: {len(result_dfs['Remains'])}")

    return result_dfs








def phmmer_random(dataframes, output_file, num_rows=10, seed=42):
    """
    Selects 10 random rows from each input Pandas DataFrame and saves them in a new DataFrame.

    :param seed: default 42, can be set manually
    :param dataframes: List of Pandas DataFrames.
    :param output_file: Path to save the output DataFrame as a CSV file.
    :param num_rows: Number of random rows to select from each DataFrame (default: 10).
    """
    selected_rows = [df.sample(n=min(num_rows, len(df)), random_state=seed) for df in dataframes]
    combined_df = pd.concat(selected_rows, ignore_index=True)
    combined_df.to_csv(output_file, index=False)
    return combined_df


def phmmer_filter_keywords(dataframe, output_prefix, keywords_dict=None):
    # Default keywords dictionary for categorization
    if not keywords_dict:
        keywords_dict = {
            "DLD": ["Dihydrolipoyl dehydrogenase",
                    "dihydrolipoyl dehydrogenase",
                    "DLD",
                    "Dihydrolipoamide dehydrogenase",
                    "dihydrolipoamide dehydrogenase"],
            "MercR": ["Mercuric reductase"],
            "GlutR": ["Glutathione reductase"]
        }

    print(f"-----------------------------------------\nSplit pHMMER hits by {keywords_dict.keys()}:")

    # Initialize a dictionary to store processed DataFrames
    results_dict = {}

    # Split results based on keywords
    split_results = split_phmmer_keywords(dataframe, keywords_dict)

    # Process each keyword and save data
    for keyword, aliases in keywords_dict.items():
        df = split_results.get(keyword)
        if df is not None:
            # Save the filtered DataFrame to a TSV file
            output_tsv_path = f"{output_prefix}{keyword.replace(' ', '_')}.tsv"
            df.to_csv(output_tsv_path, sep='\t', index=False)

            # Add the DataFrame to the results dictionary
            results_dict[keyword] = df

    # Process remaining hits
    remaining_df = split_results.get("Remains")
    if remaining_df is not None:
        # Save the remaining DataFrame to a TSV file
        remaining_tsv_path = f"{output_prefix}Remains.tsv"
        remaining_df.to_csv(remaining_tsv_path, sep='\t', index=False)

        # Add the remaining hits DataFrame to the results dictionary
        results_dict["Remains"] = remaining_df


    # Return the dictionary containing all DataFrames
    return results_dict

def phmmer_filter_E(dataframes, e_values):
    df_e = {}
    for df in dataframes:
        print(f"-----------------------------------------\nFilter pHMMER hits by E-values: {e_values}:")
        for e_value in e_values:
            filtered = filter_e_value_subset(df, e_value)
            df_e[e_value] = filtered

    return df_e


def phmmer_top_per_keyword(df_list, top=10):
    top_list = []
    for df in df_list:
        top_list.append(df.head(top))
    return top_list
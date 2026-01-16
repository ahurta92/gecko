import pandas as pd
import json
from ..dalton.daltonToJson import daltonToJson


def partition_molecule_list(mol_list):
    Flist = []
    row2 = []
    rest = []
    seconds = ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
    for mol_i in mol_list:
        try:
            if "F" in mol_i and not any([e2 in mol_i for e2 in seconds]):
                Flist.append(mol_i)
            elif any([e2 in mol_i for e2 in seconds]):
                row2.append(mol_i)
            else:
                rest.append(mol_i)
        except TypeError as f:
            print(f)
            pass

    return rest, row2, Flist


def make_detailed_df(data):
    mols = list(data.molecule.unique())
    mols = [mol for mol in mols if mol is not None]
    row1, row2, flist = partition_molecule_list(mols)
    zero = [
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "cc-pCVDZ",
        "cc-pCVTZ",
        "cc-pCVQZ",
        "cc-pV5Z",
        "cc-pV6Z",
    ]
    single = [
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
        "aug-cc-pCVDZ",
        "aug-cc-pCVTZ",
        "aug-cc-pCVQZ",
        "aug-cc-pV5Z",
        "aug-cc-pV6Z",
    ]
    double = [
        "d-aug-cc-pVDZ",
        "d-aug-cc-pVTZ",
        "d-aug-cc-pVQZ",
        "d-aug-cc-pCVDZ",
        "d-aug-cc-pCVTZ",
        "d-aug-cc-pCVQZ",
        "d-aug-cc-pV5Z",
        "d-aug-cc-pV6Z",
    ]
    triple = [
        "t-aug-cc-pVDZ",
        "t-aug-cc-pVTZ",
        "t-aug-cc-pVQZ",
        "t-aug-cc-pCVDZ",
        "t-aug-cc-pCVTZ",
        "t-aug-cc-pCVQZ",
        "t-aug-cc-pV5Z",
        "t-aug-cc-pV6Z",
    ]
    quad = [
        "q-aug-cc-pVDZ",
        "q-aug-cc-pVTZ",
        "q-aug-cc-pVQZ",
        "q-aug-cc-pCVDZ",
        "q-aug-cc-pCVTZ",
        "q-aug-cc-pCVQZ",
        "q-aug-cc-pV5Z",
        "q-aug-cc-pV6Z",
    ]

    zero_polarized = [
        "cc-pCVDZ",
        "cc-pCVTZ",
        "cc-pCVQZ",
    ]
    single_polarized = [
        "aug-cc-pCVDZ",
        "aug-cc-pCVTZ",
        "aug-cc-pCVQZ",
    ]
    double_polarized = [
        "d-aug-cc-pCVDZ",
        "d-aug-cc-pCVTZ",
        "d-aug-cc-pCVQZ",
    ]
    triple_polarized = [
        "t-aug-cc-pCVDZ",
        "t-aug-cc-pCVTZ",
        "t-aug-cc-pCVQZ",
    ]
    quadruple_polarized = [
        "q-aug-cc-pCVDZ",
        "q-aug-cc-pCVTZ",
        "q-aug-cc-pCVQZ",
    ]

    DZ = [
        "cc-pVDZ",
        "aug-cc-pVDZ",
        "d-aug-cc-pVDZ",
        "aug-cc-pCVDZ",
        "d-aug-cc-pCVDZ",
        "t-aug-cc-pVDZ",
        "q-aug-cc-pVDZ",
        "t-aug-cc-pCVDZ",
        "q-aug-cc-pCVDZ",
    ]
    TZ = [
        "cc-pVTZ",
        "aug-cc-pVTZ",
        "d-aug-cc-pVTZ",
        "aug-cc-pCVTZ",
        "d-aug-cc-pCVTZ",
        "t-aug-cc-pVTZ",
        "q-aug-cc-pVTZ",
        "t-aug-cc-pCVTZ",
        "q-aug-cc-pCVTZ",
    ]
    QZ = [
        "cc-pVQZ",
        "aug-cc-pVQZ",
        "d-aug-cc-pVQZ",
        "aug-cc-pCVQZ",
        "d-aug-cc-pCVQZ",
        "t-aug-cc-pVQZ",
        "q-aug-cc-pVQZ",
        "t-aug-cc-pCVQZ",
        "q-aug-cc-pCVQZ",
    ]
    FZ = ["cc-pV5Z", "aug-cc-pV5Z", "d-aug-cc-pV5Z"]
    SZ = ["cc-pV6Z", "aug-cc-pV6Z", "d-aug-cc-pV6Z"]

    data = data.copy()
    data["augmentation"] = ""
    data.loc[data["basis"].isin(double), "augmentation"] = "d-aug"
    data.loc[data["basis"].isin(single), "augmentation"] = "aug"
    data.loc[data["basis"].isin(triple), "augmentation"] = "t-aug"
    data.loc[data["basis"].isin(quad), "augmentation"] = "q-aug"
    data["augmentation"] = data["augmentation"].astype("category")

    data["polarization"] = "V"
    data.loc[
        data["basis"].isin(
            zero_polarized
            + single_polarized
            + double_polarized
            + triple_polarized
            + quadruple_polarized
        ),
        "polarization",
    ] = "CV"
    data["polarization"] = data["polarization"].astype("category")
    data["mol_system"] = "First-row"
    data.loc[data["molecule"].isin(row2), "mol_system"] = "Second-row"
    data["mol_system"] = data["mol_system"].astype("category")
    data["valence"] = "D"
    data.loc[data["basis"].isin(DZ), "valence"] = "D"
    data.loc[data["basis"].isin(TZ), "valence"] = "T"
    data.loc[data["basis"].isin(QZ), "valence"] = "Q"
    data.loc[data["basis"].isin(FZ), "valence"] = "5"
    data.loc[data["basis"].isin(SZ), "valence"] = "6"

    data["valence"] = data["valence"].astype("category")
    try:
        valence = list(data["valence"].unique())
        print(valence)
        data["valence"] = data["valence"].cat.reorder_categories(valence)
    except ValueError:
        data["valence"] = data["valence"].cat.reorder_categories(
            [
                "D",
                "T",
                "Q",
            ]
        )

    try:
        data["polarization"].cat.reorder_categories(
            [
                "V",
                "CV",
            ]
        )
    except ValueError as e:
        print(e)

    data["Type"] = data[["augmentation", "polarization"]].apply(
        lambda x: "-cc-p".join(x) + "nZ", axis=1
    )
    data["Type"] = data["Type"].astype("category")

    # rename type  -cc-pVnZ and -cc-pCVnZ to cc-pVnZ and cc-pCVnZ
    data["Type"] = data["Type"].cat.rename_categories(
        {
            "-cc-pVnZ": "cc-pVnZ",
            "-cc-pCVnZ": "cc-pCVnZ",
        }
    )

    possible_types = [
        "cc-pVnZ",
        "cc-pCVnZ",
        "aug-cc-pVnZ",
        "aug-cc-pCVnZ",
        "d-aug-cc-pVnZ",
        "d-aug-cc-pCVnZ",
        "t-aug-cc-pVnZ",
        "t-aug-cc-pCVnZ",
        "q-aug-cc-pVnZ",
        "q-aug-cc-pCVnZ",
    ]
    # drop if not in the actual types
    actual_types = list(data["Type"].unique())
    # drop from possible types if not in actual types
    possible_types = [t for t in possible_types if t in actual_types]
    # now reorder the categories to match the actual types
    data["Type"] = data["Type"].cat.reorder_categories(
        possible_types,
    )
    try:
        datacopy = data.copy()
        data.Type = data.Type.cat.add_categories("MRA")
        data.loc[data["basis"] == "mra-high", "Type"] = "MRA"
    except ValueError as e:
        print(e)
        return datacopy

    return data


def read_quad_basis_data(output_file, output_json, mol, basis):

    with open(output_file, "r") as daltonOutput:
        dj = daltonToJson()
        dalton_json = json.loads(dj.convert(daltonOutput))
        # get the quad data df
        quad_data = dj.readQuadResponse(output_file)
        alpha_data = dj.readAlphaFromQUADResponse(output_file)
        quad_data["molecule"] = mol
        quad_data["basis"] = basis

        alpha_data["molecule"] = mol
        alpha_data["basis"] = basis

        dalton_json["alpha"] = alpha_data.to_dict()
        dalton_json["Quad"] = quad_data.to_dict()
        # write the dalton_json to file
        with open(output_json, "w") as f:
            json.dump(dalton_json, f, indent=4)


def reconstruct_full_beta_ijk(freq_data):
    ex = freq_data.copy()
    ex["equal_ijk"] = freq_data["Beta"].apply(
        lambda x: str(x) if not isinstance(x, float) else str(None)
    )
    ex["equal_ijk"] = ex["equal_ijk"].apply(
        lambda x: ("".join(x[5:-1].split(","))) if str(None) else str(None)
    )
    ex.set_index("ijk", inplace=True)
    for index, row in ex.iterrows():
        if row.equal_ijk != str():
            ex.at[index, "Beta"] = ex.at[row.equal_ijk, "Beta"]
    ex.reset_index(inplace=True)
    ex.drop(columns="equal_ijk", inplace=True)
    return ex


def read_basis_quad(json_file):
    basis_data = pd.DataFrame.from_dict(json_file["Quad"])
    basis_data.rename(
        columns={"A-freq": "Afreq", "B-freq": "Bfreq", "C-freq": "Cfreq"}, inplace=True
    )
    # round Afreq Bfreq Cfreq to 3 decimal places
    basis_data.Afreq = basis_data.Afreq.round(3)
    basis_data.Bfreq = basis_data.Bfreq.round(3)
    basis_data.Cfreq = basis_data.Cfreq.round(3)
    # combine A-freq B-freq C-freq to ijk
    basis_data["ijk"] = (
        basis_data["A"].astype(str)
        + basis_data["B"].astype(str)
        + basis_data["C"].astype(str)
    )
    basis_data = basis_data.drop(columns=["A", "B", "C"])
    # rename Beta Value to Beta
    basis_data.rename(columns={"Beta Value": "Beta"}, inplace=True)
    bfreqs = basis_data.Bfreq.unique()
    cfreqs = basis_data.Cfreq.unique()

    ex = basis_data.copy()
    ex["equal_ijk"] = ex["Beta"].apply(
        lambda x: str(x) if not isinstance(x, float) else str(None)
    )
    ex["equal_ijk"] = ex["equal_ijk"].apply(
        lambda x: ("".join(x[5:-1].split(","))) if str(None) else str(None)
    )
    ex.set_index(["ijk", "Bfreq", "Cfreq"], inplace=True)
    ex.sort_index(inplace=True)
    ex = ex[~ex.index.duplicated(keep="first")]
    for index, row in ex.iterrows():
        if row.equal_ijk != str():
            # create new index = to equal_ijk,Cfreq,Bfreq
            new_index = (row.equal_ijk, index[2], index[1])
            ex.loc[index, "Beta"] = ex.loc[new_index, "Beta"]
    ex.reset_index(inplace=True)
    ex.drop(columns="equal_ijk", inplace=True)
    return ex


def query_beta_data(df, omega_b, omega_c):
    om = df.Cfreq.unique().sort()

    b = om[omega_b]
    c = om[omega_c]

    return df.query("Bfreq==@b and Cfreq==@c")


def query_alpha_data(df, omega):
    om = df.omega.unique()
    b = om[omega]
    return df.query("omega==@b")


def process_beta_df(beta_df):
    beta_df.Afreq = beta_df.Afreq.round(3)
    beta_df.Bfreq = beta_df.Bfreq.round(3)
    beta_df.Cfreq = beta_df.Cfreq.round(3)
    beta_df["ijk"] = (
        beta_df["A"].astype(str) + beta_df["B"].astype(str) + beta_df["C"].astype(str)
    )
    beta_df = beta_df.drop(columns=["A", "B", "C"])
    return beta_df


class BASISData:
    def __init__(self, mol, basis_set, xc, op, data_base, reread=False):
        self.mol = mol
        self.xc = xc
        self.op = op
        self.basis_set = basis_set
        # data is stored in /dalton_root/xc/mol/op/
        # and the basis set data is in quad_{mol}-{basis_set}.json
        self.data_dir = data_base
        self.output_path = self.data_dir.joinpath(
            "dalton/{}/{}/dipole/quad_{}-{}.out".format(xc, mol, mol, self.basis_set)
        )
        self.output_json_path = self.data_dir.joinpath(
            "dalton/{}/{}/dipole/quad_{}-{}.json".format(xc, mol, mol, self.basis_set)
        )

        if not self.output_json_path.exists() or reread == True:
            read_quad_basis_data(
                self.output_path, self.output_json_path, self.mol, self.basis_set
            )
            self.output_json = json.load(open(self.output_json_path))
        else:
            self.output_json = json.load(open(self.output_json_path))

        try:
            self.alpha = pd.DataFrame(self.output_json["alpha"])
            self.beta = read_basis_quad(self.output_json)
            self.energy = self.output_json["simulation"]["calculations"][1][
                "calculationResults"
            ]["totalEnergy"]["value"]
        except:
            print("error reading data")
            self.alpha = None
            self.beta = None
            self.energy = None

        # calc_path = self.data_dir.joinpath("calc_path.json")
        # self.calc_path_json = json.load(open(calc_path))
        # self.alpha_json=pd.DataFrame(json.load(open(self.data_dir.joinpath(self.calc_path_json['alpha_json_path']))))
        # self.beta_json=process_beta_df(pd.DataFrame(json.load(open(self.data_dir.joinpath(self.calc_path_json['quadratic_json_path'])))))

    def get_beta(self, omega_b, omega_c):
        return query_beta_data(self.beta, omega_b, omega_c)

    def get_alpha(self, omega):
        return query_alpha_data(self.alpha, omega)


class MRAData:
    def __init__(self, mol, data_base):
        self.mol = mol
        self.data_dir = data_base.joinpath(mol)
        calc_path = self.data_dir.joinpath("calc_path.json")
        self.calc_path_json = json.load(open(calc_path))
        self.alpha = pd.DataFrame(
            json.load(
                open(self.data_dir.joinpath(self.calc_path_json["alpha_json_path"]))
            )
        )
        self.alpha["basis"] = "MRA"
        self.alpha["molecule"] = mol
        self.beta = process_beta_df(
            pd.DataFrame(
                json.load(
                    open(
                        self.data_dir.joinpath(
                            self.calc_path_json["quadratic_json_path"]
                        )
                    )
                )
            )
        )
        self.beta["basis"] = "MRA"
        self.beta["molecule"] = mol

        self.moldft_calc_info = json.load(
            open(
                self.data_dir.joinpath(
                    self.calc_path_json["moldft"]["outfiles"]["calc_info"]
                )
            )
        )
        self.energy = self.moldft_calc_info["return_energy"]

    def get_beta(self, omega_b, omega_c):
        return query_beta_data(self.beta, omega_b, omega_c)

    def get_alpha(self, omega):
        return query_alpha_data(self.alpha, omega)


class BasisDataFrames:
    """
    Class representing a collection of basis data frames.

    Parameters:
    - database_root (str): The root directory of the database.
    - basis_sets (list): List of basis sets.
    - molecules (list): List of molecules.
    - xc (str, optional): The exchange-correlation functional. Defaults to "hf".
    - op (str, optional): The operator. Defaults to "dipole".

    Attributes:
    - database_root (str): The root directory of the database.
    - csv_data (str): The path to the CSV data directory.
    - energy (DataFrame): Data frame for energy values.
    - alpha (DataFrame): Data frame for alpha values.
    - beta (DataFrame): Data frame for beta values.
    - xc (str): The exchange-correlation functional.
    - op (str): The operator.
    - basis_sets (list): List of basis sets.
    - molecules (list): List of molecules.

    Methods:
    - collect_data(reconstruct=False): Collects the basis data data, if reconstruct is True, the data json files are reconstructed.
    - to_csv(): Saves the data frames to CSV files in the database_root / csv_files.
    - from_csv(): Loads the data frames from CSV files in database_root / csv_files.
    """

    def __init__(self, database_root, basis_sets, molecules, xc="hf", op="dipole"):
        self.database_root = database_root
        self.csv_data = database_root.joinpath("csv_data")

        if not self.csv_data.exists():
            self.csv_data.mkdir()
        self.energy = pd.DataFrame()
        self.alpha = pd.DataFrame()
        self.beta = pd.DataFrame()
        self.xc = xc
        self.op = op
        self.basis_sets = basis_sets
        self.molecules = molecules

    def collect_data(self, reconstruct=False):
        """
        Collects the basis data.

        Parameters:
        - reconstruct (bool, optional): Whether to reconstruct the data. Defaults to False.

        Returns:
        - basis_data (dict): Dictionary containing the basis data.
        - basis_available_data (dict): Dictionary containing the available basis data.
        - basis_not_available_data (dict): Dictionary containing the unavailable basis data.
        """
        basis_data = {}
        basis_available_data = {}
        basis_not_available_data = {}
        for basis in self.basis_sets:
            basis_data[basis] = {}
            basis_available_data[basis] = []
            basis_not_available_data[basis] = []
            for mol in self.molecules:
                try:
                    basis_data[basis][mol] = BASISData(
                        mol, basis, self.xc, self.op, self.database_root, reconstruct
                    )
                    basis_available_data[basis].append(mol)
                except FileNotFoundError as e:
                    basis_not_available_data[basis].append(mol)
                except:
                    basis_data[mol] = None
        self._basis_data_dict = basis_data
        self._available_data_dict = basis_available_data
        self._not_available_data_dict = basis_not_available_data
        self._from_dict(basis_data)

        return basis_data, basis_available_data, basis_not_available_data

    def _from_dict(self, data_dict):
        """
        Converts the basis data from a dictionary to data frames.

        Parameters:
        - data_dict (dict): Dictionary containing the basis data.
        """
        self.energy = pd.concat(
            [
                pd.Series(
                    {
                        "energy": data_dict[basis][mol].energy,
                        "basis": basis,
                        "molecule": mol,
                    }
                )
                for basis in self.basis_sets
                for mol in self._available_data_dict[basis]
                if data_dict[basis][mol].energy is not None
            ],
            axis=1,
        ).T
        self.alpha = pd.concat(
            [
                data_dict[basis][mol].alpha
                for basis in self.basis_sets
                for mol in self._available_data_dict[basis]
                if data_dict[basis][mol].alpha is not None
            ]
        )
        self.beta = pd.concat(
            [
                data_dict[basis][mol].beta
                for basis in self._available_data_dict
                for mol in self._available_data_dict[basis]
                if data_dict[basis][mol].beta is not None
            ]
        )

    def to_csv(self):
        """
        Saves the data frames to CSV files.
        """
        self.energy.to_csv(self.csv_data.joinpath("basis_energy.csv"), index=False)
        self.alpha.to_csv(self.csv_data.joinpath("basis_alpha.csv"), index=False)
        self.beta.to_csv(self.csv_data.joinpath("basis_beta.csv"), index=False)

    def from_csv(self):
        """
        Loads the data frames from CSV files.
        """
        self.energy = pd.read_csv(self.csv_data.joinpath("basis_energy.csv"))
        self.alpha = pd.read_csv(self.csv_data.joinpath("basis_alpha.csv"))
        self.beta = pd.read_csv(self.csv_data.joinpath("basis_beta.csv"))

    def available_data(self):
        """
        Returns the list of available data.

        Returns:
            pd.Series: The list of available data.
        """
        return self._available_data_dict

    def not_available_data(self):
        """
        Returns the list of not available data.

        Returns:
            pd.Series: The list of not available data.
        """
        return self._not_available_data_dict


import pandas as pd


class MRADataFrames:
    """
    A class that represents the data frames for MRA analysis.

    Attributes:
        database_root (str): The root directory of the database.
        molecules (list): A list of molecules.
        xc (str): The exchange-correlation functional (default is "hf").
        op (str): The operator (default is "dipole").
        energy (pd.DataFrame): The energy data frame.
        alpha (pd.DataFrame): The alpha data frame.
        beta (pd.DataFrame): The beta data frame.
    """

    def __init__(self, database_root, molecules, xc="hf", op="dipole"):
        """
        Initializes a new instance of the MRADataFrames class.

        Args:
            database_root (str): The root directory of the database.
            molecules (list): A list of molecules.
            xc (str, optional): The exchange-correlation functional (default is "hf").
            op (str, optional): The operator (default is "dipole").
        """
        self.database_root = database_root
        self.csv_data = database_root.joinpath("csv_data")
        if not self.csv_data.exists():
            self.csv_data.mkdir()
        self.molecules = molecules
        self.xc = xc
        self.op = op

        self.energy = pd.DataFrame()
        self.alpha = pd.DataFrame()
        self.beta = pd.DataFrame()

    def collect_data(self):
        """
        Collects the data for the molecules.

        Returns:
            tuple: A tuple containing the data dictionary, available data list, and not available data list.
        """
        data = {}
        output_data = self.database_root.joinpath("output")
        available_data = []
        not_available_data = []

        for mol in self.molecules:
            try:
                data[mol] = MRAData(mol, output_data)
                available_data.append(mol)
            except FileNotFoundError as e:
                not_available_data.append(mol)
            except:
                data[mol] = None
        self._from_dict(data)
        return data, available_data, not_available_data

    def _from_dict(self, data_dict):
        """
        Converts the data dictionary to data frames.

        Args:
            data_dict (dict): The data dictionary.
        """
        self.energy = pd.concat(
            [
                pd.Series(
                    {"energy": data_dict[mol].energy, "basis": "MRA", "molecule": mol}
                )
                for mol in data_dict.keys()
            ],
            axis=1,
        ).T
        self.alpha = pd.concat([data_dict[mol].alpha for mol in data_dict.keys()])
        self.beta = pd.concat([data_dict[mol].beta for mol in data_dict.keys()])

    def to_csv(self):
        """
        Saves the data frames to CSV files.
        """
        self.energy.to_csv(self.csv_data.joinpath("mra_energy.csv"), index=False)
        self.alpha.to_csv(self.csv_data.joinpath("mra_alpha.csv"), index=False)
        self.beta.to_csv(self.csv_data.joinpath("mra_beta.csv"), index=False)

    def from_csv(self):
        """
        Loads the data frames from CSV files.
        """
        self.energy = pd.read_csv(self.csv_data.joinpath("mra_energy.csv"))
        self.alpha = pd.read_csv(self.csv_data.joinpath("mra_alpha.csv"))
        self.beta = pd.read_csv(self.csv_data.joinpath("mra_beta.csv"))

    def available_data(self):
        """
        Returns the list of available data.

        Returns:
            pd.Series: The list of available data.
        """
        return self.energy.molecule.unique()

    def not_available_data(self):
        """
        Returns the list of not available data.

        Returns:
            pd.Series: The list of not available data.
        """
        return self.molecules[~self.molecules.isin(self.available_data())]


class MRAComparedBasisDF(pd.DataFrame):
    def __init__(
        self,
        polar_data,
        index,
        values: list,
        PercentError: bool,
        mra="mra-high",
        *args,
        **kwargs,
    ):

        # Use the special_parameter to modify the DataFrame or perform additional initialization
        basis_data = polar_data.query("basis!=@mra").copy()
        basis_data = basis_data.set_index(index)

        for value in values:
            basis_data[f"{value}[{mra}]"] = polar_data.query("basis==@mra").set_index(
                index
            )[value]
            if PercentError:
                basis_data[f"{value}E"] = (
                    basis_data[value] - basis_data[f"{value}[{mra}]"]
                ) / basis_data[f"{value}[{mra}]"]
            else:
                basis_data[f"{value}E"] = (
                    basis_data[value] - basis_data[f"{value}[{mra}]"]
                )
        basis_data = basis_data.reset_index()
        # create a column of percent error in alpha
        basis_data = make_detailed_df(basis_data)

        super().__init__(basis_data, *args, **kwargs)

    def get_df(self):
        return self.basis_data.copy()

    def molecule(self, mol):
        return self.query(f"molecule=='{mol}'")

    def basis(self, basis):
        return self.query(f"basis=='{basis}'")

    def freq(self, om1, om2):
        return query_beta_data(self, om1, om2)

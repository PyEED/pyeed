import io
import json
from typing import List, Dict, Any, Optional

import httpx
from Bio.SearchIO import BlastIO
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class BlastN(AbstractTool):
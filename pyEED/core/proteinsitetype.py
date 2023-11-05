from enum import Enum


class ProteinSiteType(Enum):
    ACTIVE = "active"
    BINDING = "binding"
    METAL_BINDING = "metal"
    POST_TRANS_MODIFICATION = "post-translational modification"
    UNANNOTATED = "unannotated"

    @classmethod
    def match_string(cls, s):
        for type in cls:
            if type.value == s:
                return type.value

        # print(f"unmatched site type: {s}")
        return cls.UNANNOTATED.value

from detection.app.handler.argohandler import ArgoHandler
from detection.app.handler.argpore2handler import Argpore2Handler
from detection.app.handler.deeparghandler import DeepargHandler
from detection.app.handler.kmahandler import KMAHandler
from detection.app.handler.shortbredhandler import ShortbredHandler

handler_by_key = {
    ArgoHandler.key: ArgoHandler,
    Argpore2Handler.key: Argpore2Handler,
    DeepargHandler.key: DeepargHandler,
    KMAHandler.key: KMAHandler,
    ShortbredHandler.key: ShortbredHandler,
}
